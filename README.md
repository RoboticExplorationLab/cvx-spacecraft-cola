# Convex Maneuver Planning for Spacecraft Collision Avoidance

## Table of Contents
- [About](#about)
- [Organization](#organization)
- [Dependencies](#dependecies)

## About
Conjunction analysis and maneuver planning for spacecraft collision avoidance remains a manual and time-consuming process, typically involving repeated forward simulations of hand-designed maneuvers. With the growing density of satellites in low-Earth orbit (LEO), autonomy is becoming essential for efficiently evaluating and mitigating collisions. 
In this work, we present an algorithm to design low-thrust collision-avoidance maneuvers for short-term conjunction events. We first formulate the problem as a nonconvex quadratically-constrained quadratic program (QCQP), which we then relax into a convex semidefinite program (SDP) using Shor's relaxation. We demonstrate empirically that the relaxation is tight,  which enables the recovery of globally optimal solutions to the original nonconvex problem. Our formulation produces a minimum-energy solution while ensuring a desired probability of collision at the time of closest approach. We validate our algorithm with a high-fidelity simulation of a satellite conjunction in low-Earth orbit with a simulated conjunction data message (CDM), demonstrating its effectiveness in reducing collision risk. 

## Organization
This repository contains its own Julia 1.11.6 environment specified by the Project.toml and Manifest.toml files. 

The directories in the project are the following: 
- src: algorithm source code
- examples: contains case studies from the paper

## Dependencies
The solver used for this work is Mosek. A free academic license can be obtained here: https://www.mosek.com/products/academic-licenses/