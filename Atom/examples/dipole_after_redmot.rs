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
// use lib::laser_cooling::photons_scattered::ActualPhotonsScatteredVector;
// use lib::laser_cooling::CoolingLight;
use nalgebra::Vector3;
use specs::prelude::*;
use specs::{Builder, World};
use std::time::Instant;
use rand_distr::{Normal, Distribution};
use lib::shapes::Cuboid;
use lib::sim_region::{SimulationVolume, VolumeType};
use lib::destructor::ToBeDestroyed;
use crate::laser::index::LaserIndex;
use std::f64::consts::E;
use rand_distr::num_traits::Pow;
use std::fmt;
use std::io::Write;
use std::fs;

#[derive(Debug)]
struct Ramp {
    step: i32,
    time: f64,
    rate: f64,
    delta: f64,
    initial_power: f64,
    file: String,
}

impl Component for Ramp {
    type Storage = VecStorage<Self>;
}

impl Ramp {
    pub fn update_step(&mut self, step: i32, timestep: f64) {
        self.time = (step as f64) * timestep;
        self.step = self.step + 1;
        // println!("{}", self.step)
    }

    pub fn output(&mut self, power: f64
    ) -> std::io::Result<()> {
        let mut f = fs::OpenOptions::new()
              .write(true)
              .append(true)
              .open(&self.file)
              .unwrap();
        write!(f, "{}\n", power);
        Ok(())
    }
}

impl<'a> System<'a> for Ramp {
    type SystemData = (WriteStorage<'a, GaussianBeam>);
    fn run(&mut self, (mut gaussian_beam): Self::SystemData) {
        let base: f64 = 2.0;
        let mut power = 0.0;
        for beam in (&mut gaussian_beam).join() {
            // println!("{}", self.time);
            beam.power = self.initial_power * base.pow(-1.0 * self.rate*self.time);
            power = beam.power;
        }
    Ramp::update_step(self, self.step, self.delta);
    Ramp::output(self, power);
    println!("{}", power);
    }
}

fn main() {

    let power_file_location: String = "/Users/brianbostwick/Mirror/Cambridge/QMBP_Lab/DipoleTrap_EvapCooling/EvapData/power_distro02.txt".to_string();
    let mut f = fs::File::create(&power_file_location).expect("Unable to create file");;

    let now = Instant::now();

    // Create the simulation world and builder for the ECS dispatcher.
    let mut world = World::new();
    ecs::register_components(&mut world);
    ecs::register_resources(&mut world);
    world.register::<Ramp>();
    let mut builder = ecs::create_simulation_dispatcher_builder();

    // Configure simulation output.
    builder = builder.with(
        file::new::<atom::Position, Text>("/Users/brianbostwick/Mirror/Cambridge/QMBP_Lab/DipoleTrap_EvapCooling/EvapData/pos_distro02.txt".to_string(), 10),
        "",
        &[],
    );
    builder = builder.with(
        file::new::<atom::Velocity, Text>("/Users/brianbostwick/Mirror/Cambridge/QMBP_Lab/DipoleTrap_EvapCooling/EvapData/vel_distro02.txt".to_string(), 10),
        "",
        &[],
    );

    //Beam paramters
    let power = 5.0; //initial power of the ramp
    let e_radius = 50.0e-6 / (2.0_f64.sqrt());
    let wavelength = 1064.0e-9;

    //time step
    let dt = 1.0e-6;

    //dispatcher
    let mut dispatcher = builder
        .with( Ramp { step: 0, time: 0.0, rate: 4.0_f64/3.0_f64, delta: dt, initial_power: power, file: power_file_location}, "ramp", &[])
        .build();
    dispatcher.setup(&mut world);

    //creating atom pos distrobution
    let pos_normal = Normal::new(0.0, 1.0e-6).unwrap();
    pos_normal.sample(&mut rand::thread_rng());

    //creating atom vel distrobution
    let vel_normal = Normal::new(0.0, 1.0e-3).unwrap();
    vel_normal.sample(&mut rand::thread_rng());

    for _i in 0..1000{
        world
            .create_entity()
            .with(Position {
                pos: Vector3::new((pos_normal.sample(&mut rand::thread_rng()) as f64),(pos_normal.sample(&mut rand::thread_rng()) as f64),(pos_normal.sample(&mut rand::thread_rng()) as f64)),
            })
            .with(Velocity {
                vel: Vector3::new(20.0*(vel_normal.sample(&mut rand::thread_rng()) as f64),20.0*(vel_normal.sample(&mut rand::thread_rng()) as f64),20.0*(vel_normal.sample(&mut rand::thread_rng()) as f64)),
            })
            .with(dipole::Polarizability::calculate_for(
                     wavelength, 461e-9, 2.0*3.1415*30.5e6,
                 ))
            .with(Atom)
            .with(Force::new())
            .with(AtomicTransition::strontium())
            .with(Mass{ value: 87.0 })
            .with(NewlyCreated)
            .build();
    }

    // Use a simulation bound so that atoms that escape the capture region are deleted from the simulation.
    world
        .create_entity()
        .with(Position {
            pos: Vector3::new(0.0, 0.0, 0.0),
        })
        .with(Cuboid {
            half_width: Vector3::new(5.0e-6, 5.0e-6, 5.0e-6),
        })
        .with(SimulationVolume {
            volume_type: VolumeType::Inclusive,
        })
        .build();

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

    // Define timestep
    world.insert(Timestep { delta: dt });

    // Run the simulation for a number of steps.
    for _i in 0..10000{
        dispatcher.dispatch(&mut world);
        world.maintain();
    }

    println!("Simulation completed in {} ms.", now.elapsed().as_millis());
}
