extern crate atomecs as lib;
extern crate nalgebra;
use lib::atom;
use lib::dipole;
use lib::ecs;
use lib::integrator::Timestep;
use lib::laser;
use lib::laser::gaussian::GaussianBeam;
use lib::output::file;
use lib::output::file::{Text};
use lib::atom::{Atom, AtomicTransition, Force, Mass, Position, Velocity};
use lib::atom_sources::central_creator::{CentralCreator};
use lib::atom_sources::mass::{MassDistribution};
use lib::atom_sources::emit::{AtomNumberToEmit};
use lib::atom_sources;
use lib::initiate::NewlyCreated;
use lib::laser_cooling::photons_scattered::ActualPhotonsScatteredVector;
use lib::laser_cooling::CoolingLight;
use nalgebra::Vector3;
use specs::prelude::*;
use specs::{Builder, World};
use std::time::Instant;
use rand_distr::{Normal, Distribution};
use lib::shapes::Cuboid;
use lib::sim_region::{SimulationVolume, VolumeType};

use std::fs::File;
use std::fs;
use std::io::prelude::*;
use std::io::{Write, BufReader, BufRead, Error};
use std::fs::OpenOptions;

use std::fmt;

// #[derive(Debug, Clone)]
// struct Ramp {
//     ramp: f32,
// }
//
// impl Component for Ramp {
//     type Storage = VecStorage<Self>;
// }

#[derive(Debug, Clone)]
pub struct UpdateRamp{
    ramp: f64,
}

impl Component for UpdateRamp {
	type Storage = VecStorage<Self>;
}

impl<'a> System<'a> for UpdateRamp {
    type SystemData = (
        WriteStorage<'a, GaussianBeam>
    );
    fn run(
        &mut self, ( mut power ): Self::SystemData) {
        for power in ( &mut power ).join() {

            let mut data = power.power.to_string();
            data.push_str(&"\n");
            let mut f = OpenOptions::new()
                .append(true)
                .create(true) // Optionally create the file if it doesn't already exist
                .open("power_ramp.txt")
                .expect("Unable to open file");
            f.write_all(data.as_bytes()).expect("Unable to write data");

            // if power.power > 0.0 {
            //     self.ramp = power.power - 0.0008;
            //     power.power = self.ramp;
            // }
            // else {
            //     power.power = 0.0;
            // }
            println!("Power: {:?}", i);
        }
    }
}


// impl fmt::Display for UpdateRamp {
// 	fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
// 		write!(f, "{:?}", self.ramp)
// 	}
// }

fn main() {

    let now = Instant::now();

    // Create the simulation world and builder for the ECS dispatcher.
    let mut world = World::new();
    ecs::register_components(&mut world);
    ecs::register_resources(&mut world);
    world.register::<UpdateRamp>();
    let mut builder = ecs::create_simulation_dispatcher_builder();

    // Configure simulation output.
    builder = builder.with(
        file::new::<atom::Position, Text>("pos_dipole.txt".to_string(), 10),
        "",
        &[],
    );
    builder = builder.with(
        file::new::<atom::Velocity, Text>("vel_dipole.txt".to_string(), 10
    ),
        "",
        &[],
    );
    // builder = builder.with(
    //     file::new::<UpdateRamp, Text>("power.txt".to_string(), 10),
    //     "",
    //     &[],
    // );

    //Beam paramters
    let power = 8.0; //initial power of the ramp
    let e_radius = 60.0e-6 / (2.0_f64.sqrt());
    let wavelength = 1064.0e-9;
    let val = power;


    //dispatcher
    let mut dispatcher = builder
        .with(UpdateRamp { ramp: val }, "ramp", &[])
        .build();
    dispatcher.setup(&mut world);

    //creating atom pos distrobution
    let pos_normal = Normal::new(0.0, 6000.0e-6 / (2.0_f64.sqrt())).unwrap();
    pos_normal.sample(&mut rand::thread_rng());

    //creating atom vel distrobution
    let vel_normal = Normal::new(0.0, 10.0).unwrap();
    vel_normal.sample(&mut rand::thread_rng());

    for i in 0..100 {
        world
            .create_entity()
            .with(Position {
                pos: Vector3::new((pos_normal.sample(&mut rand::thread_rng()) as f64)*10.0e-6,(pos_normal.sample(&mut rand::thread_rng()) as f64)*1.0e-6,(pos_normal.sample(&mut rand::thread_rng()) as f64)*1.0e-6),
            })
            .with(Velocity {
                vel: Vector3::new(0.0,0.0,0.0),
            })
            .with(dipole::Polarizability::calculate_for(
                     wavelength, 461e-9, 32.0e6,
                 ))
            .with(Atom)
            .with(NewlyCreated)
            .with(Force::new())
            .with(AtomicTransition::strontium())
            .with(Mass{ value: 87.0 })
            .build();
    }

    // Create dipole laser.
    let mut gaussian_beam = GaussianBeam {
        intersection: Vector3::new(0.0, 0.0, 0.0),
        e_radius: e_radius,
        power: power,
        direction: Vector3::y(),
        rayleigh_range: crate::laser::gaussian::calculate_rayleigh_range(&wavelength, &e_radius),
        ellipticity: 0.0,
    };
    world
        .create_entity()
        .with(gaussian_beam)
        // xs.with(Ramp {ramp: 2.0} )
        .with(dipole::DipoleLight {
            wavelength: wavelength,
        })
        .with(laser::frame::Frame {
            x_vector: Vector3::x(),
            y_vector: Vector3::z(),
        })
        .build();

    let mut gaussian_beam = GaussianBeam {
        intersection: Vector3::new(0.0, 0.0, 0.0),
        e_radius: e_radius,
        power: power,
        direction: Vector3::x(),
        rayleigh_range: crate::laser::gaussian::calculate_rayleigh_range(&wavelength, &e_radius),
        ellipticity: 0.0,
    };
    world
        .create_entity()
        .with(gaussian_beam)
        .with(dipole::DipoleLight {
            wavelength: wavelength,
        })
        .with(laser::frame::Frame {
            x_vector: Vector3::y(),
            y_vector: Vector3::z(),
        })
        .build();


    // Use a simulation bound so that atoms that escape the capture region are deleted from the simulation.
    world
        .create_entity()
        .with(Position {
            pos: Vector3::new(0.0, 0.0, 0.0),
        })
        .with(Cuboid {
            half_width: Vector3::new(e_radius/4.0, e_radius/4.0, e_radius/4.0),
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();




    // Define timestep
    world.insert(Timestep { delta: 1.0e-5 });

    // Run the simulation for a number of steps.
    for _i in 0..10000{
        dispatcher.dispatch(&mut world);
        world.maintain();
    }

    println!("Simulation completed in {} ms.", now.elapsed().as_millis());
}
