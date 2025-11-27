# Time-Aware Probabilistic Modeling of IEEE 3-Bus Protection System

This repository contains a PRISM MDP model of an IEEE 3-bus distribution test system equipped with six digital overcurrent relays. The model integrates IEC normal inverse-time curves, global fault-age timing, coordination time margin (CTM) analysis, communication failures, and an auxiliary supervisory service.

The model is designed to support **time-aware quantitative verification** of protection performance, including:

- Probability of successful line isolation  
- Probability of isolation failure
- Probability of system under risk
- Probability of false tripping  
- Probability that the system remains under fault past the thermal limit  

## Files

- `ieee3bus_timed.prism`  
  Main PRISM MDP model of the protection system:
  - Fault process (`FAULTS` module)
  - Global clock (`CLOCK` module) using 50 ms time steps
  - IEC inverse-time relay trip computations
  - Six relay modules (`Relay_R1` as template + renamings for R2–R6)
  - CTM modules (`CTM_R1`–`CTM_R6`)
  - Signal dispatching modules
  - Communication channel modules between relays and breakers
  - Supervisory service (`Sup_serv1`)
  - Rewards and labels for property specification

## Requirements

- [PRISM model checker](https://www.prismmodelchecker.org/) with support for:
  - MDPs
  - PCTL properties
  - Rewards

## Example Properties

You can either embed properties in PRISM or keep them in a separate `.props` file. Examples:

```prism
// Probability of eventual isolation of each line
Pmax=? [ F "L1_Isolated" ]
Pmax=? [ F "L2_Isolated" ]
Pmax=? [ F "L3_Isolated" ]

// Probability of isolation failure
Pmax=? [ F "L1_Isol_Failed" ]
Pmax=? [ F "L2_Isol_Failed" ]
Pmax=? [ F "L3_Isol_Failed" ]

// Probability of Risk
Pmax=? [F("L1_Risk")]
Pmax=? [F("L2_Risk")]
Pmax=? [F("L3_Risk")]

// Probability of False tripping
Pmax=? [F("L1_False_Trip")]
Pmax=? [F("L2_False_Trip")]
Pmax=? [F("L3_False_Trip")] 

// Probability that any line remains under fault past thermal limit
Pmax=? [ F "FaultPastThermal" ]

// How to Run
# Example: compute max probability that L1 is eventually isolated
prism ieee3bus_timed.prism -prop 'Pmax=? [ F "L1_Isolated" ]'

# Or using an external properties file:
prism ieee3bus_timed.prism props.props -prop 1



