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
use std::f64::consts::PI;
use lib::integrator::INTEGRATE_VELOCITY_SYSTEM_NAME;

use atomecs::output::file::SerdeJson;

fn main() {

    let now = Instant::now();

    // Create the simulation world and builder for the ECS dispatcher.
    let mut world = World::new();
    ecs::register_components(&mut world);
    ecs::register_resources(&mut world);
    let mut builder = ecs::create_simulation_dispatcher_builder();

    let mut loc1 = "/Users/brianbostwick/Mirror/Cambridge/QMBP_Lab/DipoleTrap_EvapCooling/data/range.txt";
    // let mut loc2 = "/Users/brianbostwick/Mirror/Cambridge/QMBP_Lab/DipoleTrap_EvapCooling/data/range.txt";
    // let mut loc3 = "/Users/brianbostwick/Mirror/Cambridge/QMBP_Lab/DipoleTrap_EvapCooling/data/range.txt";

    // // Configure simulation output.
    // builder = builder.with(
    //     file::new::<atom::Position, Text>( loc1.to_string(), 10),
    //     "",
    //     &[],
    // );
    // builder = builder.with(
    //     file::new::<atom::Velocity, Text>( loc2.to_string(), 10),
    //     "",
    //     &[],
    // );

    builder = builder.with(
        file::new::<gaussian::get_gaussian_beam_intensity, Text>( loc1.to_string(), 10),
        "",
        &[],
    );

    //Beam paramters
    let power = 10.0; //initial power of the ramp
    let e_radius = 50.0e-6 / (2.0_f64.sqrt());
    let wavelength = 1064.0e-9;

    //dispatcher
    let mut dispatcher = builder
        .build();
    dispatcher.setup(&mut world);

    world
        .create_entity()
        .with(Position {
            pos: Vector3::new(1.0e-6, 1.0e-6, 0.0),
        })
        .with(atom::Velocity {
            vel: Vector3::new(0.0, 0.0, 0.0),
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

    // Define timestep
    world.insert(Timestep { delta: 1.0e-6 });

    // Run the simulation for a number of steps.
    for _i in 0..4000{
        dispatcher.dispatch(&mut world);
        world.maintain();
    }

    println!("Simulation completed in {} ms.", now.elapsed().as_millis());
}
