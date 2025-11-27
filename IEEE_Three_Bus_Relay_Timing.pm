mdp


// ================================================================
// IEEE 3-Bus Test System with 6 Overcurrent Relays
// Time-aware MDP model with IEC NI inverse-time curves
// ================================================================

// ---------------------------------------------------------------
// Global timing & probabilities
// ---------------------------------------------------------------
const int MAX_TIME = 16; // 800 ms / 50 ms = 16 steps
const int T_WAIT   = 16; // Supervisory thermal threshold

const double flt = 0.1;  // Fault occurrence probability
const double IED = 0.1;  // Relay internal failure probability
const double COM = 0.1;  // Communication failure probability
const double WD  = 0.1;  // Watchdog / internal error probability

// ---------------------------------------------------------------
// IEC Normal Inverse Curve parameters
// ---------------------------------------------------------------
const double K = 0.14;
const double a = 0.02;

// ---------------------------------------------------------------
// Relay settings (TMS, CTR, PS, pickup currents)
// ---------------------------------------------------------------
// Relay 1
const double TMS1 = 0.1070;
const double CTR1 = 60.0;     // 300/5 CT ratio
const double PS1  = 2.5;      // Plug setting in A (secondary)
const double Ipu_R1 = PS1 * CTR1;

// Relay 2
const double TMS2 = 0.1080;
const double CTR2 = 40.0;
const double PS2  = 2.0;
const double Ipu_R2 = PS2 * CTR2;

// Relay 3
const double TMS3 = 0.1000;
const double CTR3 = 40.0;
const double PS3  = 3.0;
const double Ipu_R3 = PS3 * CTR3;

// Relay 4
const double TMS4 = 0.1000;
const double CTR4 = 60.0;
const double PS4  = 2.5;
const double Ipu_R4 = PS4 * CTR4;

// Relay 5
const double TMS5 = 0.1000;
const double CTR5 = 40.0;
const double PS5  = 2.5;
const double Ipu_R5 = PS5 * CTR5;

// Relay 6
const double TMS6 = 0.1120;
const double CTR6 = 80.0;
const double PS6  = 1.5;
const double Ipu_R6 = PS6 * CTR6;

// ---------------------------------------------------------------
// Fault currents for primary / backup relays
// ---------------------------------------------------------------
// Primaries
const double IFp_R1 = 1978.9;
const double IFp_R2 = 1525.7;
const double IFp_R3 = 1683.9;
const double IFp_R4 = 1815.4;
const double IFp_R5 = 1499.66;
const double IFp_R6 = 1766.3;

// Backups
const double IFb_R1 = 617.22;
const double IFb_R2 = 145.34;
const double IFb_R3 = 384.0;
const double IFb_R4 = 545.0;
const double IFb_R5 = 175.0;
const double IFb_R6 = 466.17;

// ---------------------------------------------------------------
// Trip time calculations (IEC NI inverse-time, in 50 ms steps)
// ---------------------------------------------------------------
const int Ttrip_pR1 = ceil((1000 * (TMS1 * K / (pow(IFp_R1 / Ipu_R1, a) - 1))) / 50);
const int Ttrip_bR1 = ceil((1000 * (TMS1 * K / (pow(IFb_R1 / Ipu_R1, a) - 1))) / 50);

const int Ttrip_pR2 = ceil((1000 * (TMS2 * K / (pow(IFp_R2 / Ipu_R2, a) - 1))) / 50);
const int Ttrip_bR2 = ceil((1000 * (TMS2 * K / (pow(IFb_R2 / Ipu_R2, a) - 1))) / 50);

const int Ttrip_pR3 = ceil((1000 * (TMS3 * K / (pow(IFp_R3 / Ipu_R3, a) - 1))) / 50);
const int Ttrip_bR3 = ceil((1000 * (TMS3 * K / (pow(IFb_R3 / Ipu_R3, a) - 1))) / 50);

const int Ttrip_pR4 = ceil((1000 * (TMS4 * K / (pow(IFp_R4 / Ipu_R4, a) - 1))) / 50);
const int Ttrip_bR4 = ceil((1000 * (TMS4 * K / (pow(IFb_R4 / Ipu_R4, a) - 1))) / 50);

const int Ttrip_pR5 = ceil((1000 * (TMS5 * K / (pow(IFp_R5 / Ipu_R5, a) - 1))) / 50);
const int Ttrip_bR5 = ceil((1000 * (TMS5 * K / (pow(IFb_R5 / Ipu_R5, a) - 1))) / 50);

const int Ttrip_pR6 = ceil((1000 * (TMS6 * K / (pow(IFp_R6 / Ipu_R6, a) - 1))) / 50);
const int Ttrip_bR6 = ceil((1000 * (TMS6 * K / (pow(IFb_R6 / Ipu_R6, a) - 1))) / 50);

// ---------------------------------------------------------------
// CTM (Coordination Time Margin) calculations
// ---------------------------------------------------------------
formula CTM_valR1 = Ttrip_bR5 - Ttrip_pR1;
formula CTM_valR2 = Ttrip_bR4 - Ttrip_pR2;
formula CTM_valR3 = Ttrip_bR1 - Ttrip_pR3;
formula CTM_valR4 = Ttrip_bR6 - Ttrip_pR4;
formula CTM_valR5 = Ttrip_bR3 - Ttrip_pR5;
formula CTM_valR6 = Ttrip_bR2 - Ttrip_pR6;

// 1 = OK, 2 = too small, 3 = too large (300–400 ms OK)
const int CTM_R1 = (CTM_valR1 >= 6 & CTM_valR1 <= 8) ? 1 : (CTM_valR1 < 6 ? 2 : 3);
const int CTM_R2 = (CTM_valR2 >= 6 & CTM_valR2 <= 8) ? 1 : (CTM_valR2 < 6 ? 2 : 3);
const int CTM_R3 = (CTM_valR3 >= 6 & CTM_valR3 <= 8) ? 1 : (CTM_valR3 < 6 ? 2 : 3);
const int CTM_R4 = (CTM_valR4 >= 6 & CTM_valR4 <= 8) ? 1 : (CTM_valR4 < 6 ? 2 : 3);
const int CTM_R5 = (CTM_valR5 >= 6 & CTM_valR5 <= 8) ? 1 : (CTM_valR5 < 6 ? 2 : 3);
const int CTM_R6 = (CTM_valR6 >= 6 & CTM_valR6 <= 8) ? 1 : (CTM_valR6 < 6 ? 2 : 3);

// ================================================================
// Line isolation / failure formulas
// ================================================================
formula L1_Isol =
  FC1 = 1 &
  (((Break1 = true & cb1 = 1) | (Break5 = true & cb5 = 1)) &
   ((Break2 = true & cb2 = 1) | (Break4 = true & cb4 = 1)) &
   (!(sv5 = 0 & cb5 = 2) & !(sv4 = 0 & cb4 = 2)));

formula L2_Isol =
  FC2 = 1 &
  (((Break3 = true & cb3 = 1) | (Break1 = true & cb1 = 1)) &
   ((Break4 = true & cb4 = 1) | (Break6 = true & cb6 = 1)) &
   (!(sv1 = 0 & cb1 = 2) & !(sv6 = 0 & cb6 = 2)));

formula L3_Isol =
  FC3 = 1 &
  (((Break5 = true & cb5 = 1) | (Break3 = true & cb3 = 1)) &
   ((Break6 = true & cb6 = 1) | (Break2 = true & cb2 = 1)) &
   (!(sv3 = 0 & cb3 = 2) & !(sv2 = 0 & cb2 = 2)));

formula L1_Fail =
  FC1 = 1 &
  (((sv5 = 0 & WD1 = 1 & WD5 = 2 & cb5 = 2) |
    (sv4 = 0 & WD2 = 1 & WD4 = 2 & cb4 = 2)) |
   ((sv5 = 0 & cb5 = 2) | (sv4 = 0 & cb4 = 2)));

formula L2_Fail =
  FC2 = 1 &
  (((sv1 = 0 & WD3 = 1 & WD1 = 2 & cb1 = 2) |
    (sv6 = 0 & WD4 = 1 & WD6 = 2 & cb6 = 2)) |
   ((sv1 = 0 & cb1 = 2) | (sv6 = 0 & cb6 = 2)));

formula L3_Fail =
  FC3 = 1 &
  (((sv3 = 0 & WD5 = 1 & WD3 = 2 & cb3 = 2) |
    (sv2 = 0 & WD6 = 1 & WD2 = 2 & cb2 = 2)) |
   ((sv3 = 0 & cb3 = 2) | (sv2 = 0 & cb2 = 2)));

// ================================================================
// Supervisory conditions (per relay)
// ================================================================
formula sup_R1_cond1 = FC2 = 1 & (R3 = 2 | cb3 = 2) & Request3 = true & c31 = 2 & R1 = 3;
formula sup_R1_cond2 = FC2 = 1 & (R3 = 2 | cb3 = 2) & (R1 = 2);
formula sup_R1_cond3 = FC2 = 1 & (R1 = 4 | R1 = 5) & Request3 = true & c31 = 2;

formula sup_R2_cond1 = FC3 = 1 & (R6 = 2 | cb6 = 2) & Request6 = true & c62 = 2 & R2 = 3;
formula sup_R2_cond2 = FC3 = 1 & (R6 = 2 | cb6 = 2) & (R2 = 2);
formula sup_R2_cond3 = FC3 = 1 & (R2 = 4 | R2 = 5) & Request6 = true & c62 = 2;

formula sup_R3_cond1 = FC3 = 1 & (R5 = 2 | cb5 = 2) & Request5 = true & c53 = 2 & R3 = 3;
formula sup_R3_cond2 = FC3 = 1 & (R5 = 2 | cb5 = 2) & (R3 = 2);
formula sup_R3_cond3 = FC3 = 1 & (R3 = 4 | R3 = 5) & Request5 = true & c53 = 2;

formula sup_R4_cond1 = FC1 = 1 & (R2 = 2 | cb2 = 2) & Request2 = true & c24 = 2 & R4 = 3;
formula sup_R4_cond2 = FC1 = 1 & (R2 = 2 | cb2 = 2) & (R4 = 2);
formula sup_R4_cond3 = FC1 = 1 & (R4 = 4 | R4 = 5) & Request2 = true & c24 = 2;

formula sup_R5_cond1 = FC1 = 1 & (R1 = 2 | cb1 = 2) & Request1 = true & c15 = 2 & R5 = 3;
formula sup_R5_cond2 = FC1 = 1 & (R1 = 2 | cb1 = 2) & (R5 = 2);
formula sup_R5_cond3 = FC1 = 1 & (R5 = 4 | R5 = 5) & Request1 = true & c15 = 2;

formula sup_R6_cond1 = FC2 = 1 & (R4 = 2 | cb4 = 2) & Request4 = true & c46 = 2 & R6 = 3;
formula sup_R6_cond2 = FC2 = 1 & (R4 = 2 | cb4 = 2) & (R6 = 2);
formula sup_R6_cond3 = FC2 = 1 & (R6 = 4 | R6 = 5) & Request4 = true & c46 = 2;

// ================================================================
// Timing core:
//   - Primaries use tglobl (absolute fault age).
//   - Backups also keyed to tglobl (scenario-specific currents).
// ================================================================

// Max primary trip time across all relays
formula Max_p =
  max(Ttrip_pR1,
  max(Ttrip_pR2,
  max(Ttrip_pR3,
  max(Ttrip_pR4,
  max(Ttrip_pR5, Ttrip_pR6)))));

// Max backup trip time across all relays
formula Max_b =
  max(Ttrip_bR1,
  max(Ttrip_bR2,
  max(Ttrip_bR3,
  max(Ttrip_bR4,
  max(Ttrip_bR5, Ttrip_bR6)))));

// Overall upper bound for timers
formula Max_Time = max(MAX_TIME, max(Max_p, Max_b));

// Global backup timer variables (currently unused but reserved)
global Timer_bR1 : [0..Max_Time] init 0;
global Timer_bR2 : [0..Max_Time] init 0;
global Timer_bR3 : [0..Max_Time] init 0;
global Timer_bR4 : [0..Max_Time] init 0;
global Timer_bR5 : [0..Max_Time] init 0;
global Timer_bR6 : [0..Max_Time] init 0;

// ---------------------------------------------------------------
// Trip-due definitions (primary, using tglobl)
// ---------------------------------------------------------------
formula TripDue_pR1 = Arm_pR1 & (R1 = 5) & FC1 = 1 & (tglobl >= Ttrip_pR1);
formula TripDue_pR2 = Arm_pR2 & (R2 = 5) & FC1 = 1 & (tglobl >= Ttrip_pR2);
formula TripDue_pR3 = Arm_pR3 & (R3 = 5) & FC2 = 1 & (tglobl >= Ttrip_pR3);
formula TripDue_pR4 = Arm_pR4 & (R4 = 5) & FC2 = 1 & (tglobl >= Ttrip_pR4);
formula TripDue_pR5 = Arm_pR5 & (R5 = 5) & FC3 = 1 & (tglobl >= Ttrip_pR5);
formula TripDue_pR6 = Arm_pR6 & (R6 = 5) & FC3 = 1 & (tglobl >= Ttrip_pR6);

// ---------------------------------------------------------------
// Trip-due definitions (backup, also using tglobl)
// ---------------------------------------------------------------
formula TripDue_bR1 = Arm_bR1 & (R1 = 4 | R1 = 5) & (tglobl >= Ttrip_bR1);
formula TripDue_bR2 = Arm_bR2 & (R2 = 4 | R2 = 5) & (tglobl >= Ttrip_bR2);
formula TripDue_bR3 = Arm_bR3 & (R3 = 4 | R3 = 5) & (tglobl >= Ttrip_bR3);
formula TripDue_bR4 = Arm_bR4 & (R4 = 4 | R4 = 5) & (tglobl >= Ttrip_bR4);
formula TripDue_bR5 = Arm_bR5 & (R5 = 4 | R5 = 5) & (tglobl >= Ttrip_bR5);
formula TripDue_bR6 = Arm_bR6 & (R6 = 4 | R6 = 5) & (tglobl >= Ttrip_bR6);

// Any trip due across all relays
formula ANY_TRIP_DUE_RAW =
  TripDue_pR1 | TripDue_bR1 |
  TripDue_pR2 | TripDue_bR2 |
  TripDue_pR3 | TripDue_bR3 |
  TripDue_pR4 | TripDue_bR4 |
  TripDue_pR5 | TripDue_bR5 |
  TripDue_pR6 | TripDue_bR6;

// Trips are only considered "due" before thermal limit
formula ANY_TRIP_DUE =
  (tglobl < T_WAIT) & ANY_TRIP_DUE_RAW;

// Any primary/backup timer still "waiting" (based on tglobl vs trip times)
formula ANY_TIMER_ACTIVE =
     (Arm_pR1 & tglobl < Ttrip_pR1) |
     (Arm_pR2 & tglobl < Ttrip_pR2) |
     (Arm_pR3 & tglobl < Ttrip_pR3) |
     (Arm_pR4 & tglobl < Ttrip_pR4) |
     (Arm_pR5 & tglobl < Ttrip_pR5) |
     (Arm_pR6 & tglobl < Ttrip_pR6) |
     (Arm_bR1 & tglobl < Ttrip_bR1) |
     (Arm_bR2 & tglobl < Ttrip_bR2) |
     (Arm_bR3 & tglobl < Ttrip_bR3) |
     (Arm_bR4 & tglobl < Ttrip_bR4) |
     (Arm_bR5 & tglobl < Ttrip_bR5) |
     (Arm_bR6 & tglobl < Ttrip_bR6);

// ================================================================
// Supervisory due / condition formulas
// ================================================================
formula SUP_DUE1 =
  sv1 = 0 & FC2 = 1 & (sup_R1_cond1 | sup_R1_cond2 | sup_R1_cond3)
  & (tglobl >= T_WAIT);

formula SUP_DUE2 =
  sv2 = 0 & FC3 = 1 & (sup_R2_cond1 | sup_R2_cond2 | sup_R2_cond3)
  & (tglobl >= T_WAIT);

formula SUP_DUE3 =
  sv3 = 0 & FC3 = 1 & (sup_R3_cond1 | sup_R3_cond2 | sup_R3_cond3)
  & (tglobl >= T_WAIT);

formula SUP_DUE4 =
  sv4 = 0 & FC1 = 1 & (sup_R4_cond1 | sup_R4_cond2 | sup_R4_cond3)
  & (tglobl >= T_WAIT);

formula SUP_DUE5 =
  sv5 = 0 & FC1 = 1 & (sup_R5_cond1 | sup_R5_cond2 | sup_R5_cond3)
  & (tglobl >= T_WAIT);

formula SUP_DUE6 =
  sv6 = 0 & FC2 = 1 & (sup_R6_cond1 | sup_R6_cond2 | sup_R6_cond3)
  & (tglobl >= T_WAIT);

// Thermal supervisory: absolute fault age reached T_WAIT
formula ANY_THERMAL_SUP_DUE =
  Fault & (FC1 = 1 | FC2 = 1 | FC3 = 1) & (tglobl >= T_WAIT);

formula ANY_SUP_DUE =
  ANY_THERMAL_SUP_DUE |
  SUP_DUE1 | SUP_DUE2 | SUP_DUE3 | SUP_DUE4 | SUP_DUE5 | SUP_DUE6;

formula SUP_COND1 = sv1 = 0 & FC2 = 1 & (sup_R1_cond1 | sup_R1_cond2 | sup_R1_cond3);
formula SUP_COND2 = sv2 = 0 & FC3 = 1 & (sup_R2_cond1 | sup_R2_cond2 | sup_R2_cond3);
formula SUP_COND3 = sv3 = 0 & FC3 = 1 & (sup_R3_cond1 | sup_R3_cond2 | sup_R3_cond3);
formula SUP_COND4 = sv4 = 0 & FC1 = 1 & (sup_R4_cond1 | sup_R4_cond2 | sup_R4_cond3);
formula SUP_COND5 = sv5 = 0 & FC1 = 1 & (sup_R5_cond1 | sup_R5_cond2 | sup_R5_cond3);
formula SUP_COND6 = sv6 = 0 & FC2 = 1 & (sup_R6_cond1 | sup_R6_cond2 | sup_R6_cond3);

formula ANY_SUP_COND_DUE =
   SUP_COND1 | SUP_COND2 | SUP_COND3 | SUP_COND4 | SUP_COND5 | SUP_COND6;

// Clock ticks when a fault is present and some timer/supervisory condition is pending
formula NEED_TICK =
  Fault & (ANY_TIMER_ACTIVE | (ANY_SUP_COND_DUE & tglobl < T_WAIT));

// ---------------------------------------------------------------
// Worst-case line operation times and "bad coordination" flags
// ---------------------------------------------------------------
formula Max_L1_R1R5 = max(Ttrip_pR1, Ttrip_bR5);
formula Max_L1_R2R4 = max(Ttrip_pR2, Ttrip_bR4);

formula Max_L2_R3R1 = max(Ttrip_pR3, Ttrip_bR1);
formula Max_L2_R4R6 = max(Ttrip_pR4, Ttrip_bR6);

formula Max_L3_R5R3 = max(Ttrip_pR5, Ttrip_bR3);
formula Max_L3_R6R2 = max(Ttrip_pR6, Ttrip_bR2);

formula BAD_COORD_L1_sv5 = Max_L1_R1R5 >= T_WAIT;
formula BAD_COORD_L1_sv4 = Max_L1_R2R4 >= T_WAIT;
formula BAD_COORD_L2_sv1 = Max_L2_R3R1 >= T_WAIT;
formula BAD_COORD_L2_sv6 = Max_L2_R4R6 >= T_WAIT;
formula BAD_COORD_L3_sv3 = Max_L3_R5R3 >= T_WAIT;
formula BAD_COORD_L3_sv2 = Max_L3_R6R2 >= T_WAIT;

// ================================================================
// Fault process
// ================================================================
module FAULTS
  Fault : bool init false;
  FC1   : [0..2] init 0;
  FC2   : [0..2] init 0;
  FC3   : [0..2] init 0;

  // Fault arrival (single fault activated with probability flt)
  [] Fault = false & FC1 = 0 & FC2 = 0 & FC3 = 0 ->
       flt : (Fault' = true) + (1 - flt) : (Fault' = false);

  // Fault type selection (FC1, FC2, FC3) – here FC3 only (can be adapted)
  [] Fault = true & FC1 = 0 & FC2 = 0 & FC3 = 0 ->
       0 : (FC1' = 1) + 0 : (FC2' = 1) + 1 : (FC3' = 1);

  [FC1_clrd] FC1 = 1 & L1_Isol -> (FC1' = 2) & (Fault' = false);
  [FC2_clrd] FC2 = 1 & L2_Isol -> (FC2' = 2) & (Fault' = false);
  [FC3_clrd] FC3 = 1 & L3_Isol -> (FC3' = 2) & (Fault' = false);
endmodule

// ================================================================
// Global clock
//   - tglobl encodes absolute fault age in 50 ms steps
// ================================================================
module CLOCK
  tglobl  : [0..Max_Time] init 0;
  started : bool init false;

  [start] !started & (FC1 = 1 | FC2 = 1 | FC3 = 1) ->
    (started' = true);

  [] started & NEED_TICK
     & !ANY_TRIP_DUE & !ANY_SUP_DUE
     & tglobl < Max_Time ->
       (tglobl' = tglobl + 1);
endmodule

// ================================================================
// Relay R1 (template for all relays via renaming)
// ================================================================
// R1 is main for FC1 (backup on L3 via R5) and backup for FC2 / R3
module Relay_R1
  // R1 states: 0 Idle, 1 Trip, 2 Fail, 3 Lockout, 4 Reset, 5 Active
  R1 : [0..5] init 0;
  // Watchdog: 0 idle, 1 Error, 2 No Error
  WD1 : [0..2] init 0;
  R1_locked : bool init false;

  // Arm flags
  Arm_pR1 : bool init false;
  Arm_bR1 : bool init false;

  // Watchdog check: relay becomes active or fails
  [] R1 = 0 & WD1 = 0 & (FC1 = 1 | FC2 = 1) ->
     1 - WD : (WD1' = 2) & (R1' = 5)
   +     WD : (WD1' = 1) & (R1' = 2);

  // R1 as PRIMARY for FC1 (arming just enables/gates)
  [] FC1 = 1 & R1 = 5 & (CT1 > 0 | WD5 = 1) & Lock1 = true & !Arm_pR1 ->
     (Arm_pR1' = true);

  // R1 as BACKUP for FC2 & R3
  [] !Arm_bR1 & FC2 = 1 & (R1 = 5 | R1 = 4) & Request3 = true & c31 = 1 &
     ((WD3 = 2 & R3 = 2 & t3 = 2) |
      (cb3 = 2 & R3 = 1 & t3 = 3) |
      (WD3 = 1 & R3 = 2 & t3 = 4)) ->
     (Arm_bR1' = true);

  // R1 as PRIMARY for FC1 – trip when global time reaches Ttrip_pR1
  [Trp_prm1] Arm_pR1 & R1 = 5 & FC1 = 1 & (CT1 > 0 | WD5 = 1) & Lock1 = true
              & (tglobl >= Ttrip_pR1) ->
              1 - IED : (R1' = 1) & (Arm_pR1' = false)
            +     IED : (R1' = 2) & (Arm_pR1' = false);

  // R1 as BACKUP for FC2 & R3 – lockout, reset, and trip on request
  [Lock_r1]  R1 = 5 & (CT3 = 2 | CT3 = 3) & Lock3 = true & c31 = 1 & FC2 = 1 & t3 = 1 ->
              (R1' = 3);

  [Lock_r1]  R1 = 5 & (CT3 = 1) & Lock3 = true & c31 = 1 & FC2 = 1 & t3 = 1 ->
              (R1' = 3);

  // False trip when CT3 too small and return channel fails
  [False_trip1] R1 = 5 & CT3 = 2 & Lock3 = true & c31 = 2 & FC2 = 1 & Request3 = false
                  & (tglobl >= Ttrip_bR1) ->
                  1 - IED : (R1' = 1) & (Arm_bR1' = false)
                +     IED : (R1' = 2) & (Arm_bR1' = false);

  [Reset_r1]  R1 = 3 & (CT3 > 0) & Reset3 = true & c31 = 1 & FC2 = 1 & t3 = 2 ->
              (R1' = 4);

  [Locked_r1] R1 = 3 & (CT3 > 0) & Reset3 = true & c31 = 2 & FC2 = 1 & t3 = 2
              & R1_locked = false ->
              (R1' = 3) & (R1_locked' = true);

  // Backup trip of R1 after its backup delay, on valid request from R3
  [Trp_bkp1] (R1 = 4) & Request3 = true & c31 = 1 & FC2 = 1 & WD3 = 2 & R3 = 2 & t3 = 2
              & (tglobl >= Ttrip_bR1) ->
              1 - IED : (R1' = 1) & (Arm_bR1' = false)
            +     IED : (R1' = 2) & (Arm_bR1' = false);

  [Trp_bkp1] (R1 = 4) & Request3 = true & c31 = 1 & FC2 = 1 & cb3 = 2 & R3 = 1 & t3 = 3
              & (tglobl >= Ttrip_bR1) ->
              1 - IED : (R1' = 1) & (Arm_bR1' = false)
            +     IED : (R1' = 2) & (Arm_bR1' = false);

  [Trp_bkp1] (R1 = 5) & Request3 = true & c31 = 1 & FC2 = 1 & WD3 = 1 & R3 = 2 & t3 = 4
              & (tglobl >= Ttrip_bR1) ->
              1 - IED : (R1' = 1) & (Arm_bR1' = false)
            +     IED : (R1' = 2) & (Arm_bR1' = false);
endmodule

// ================================================================
// Relay renamings for R2–R6
// ================================================================
module Relay_R2 = Relay_R1 [
  R1=R2, WD1=WD2, R1_locked=R2_locked, WD5=WD4, Lock1=Lock2,
  Timer_bR1=Timer_bR2, Ttrip_pR1=Ttrip_pR2, Ttrip_bR1=Ttrip_bR2,
  FC1=FC1, CT1=CT2, FC2=FC3, CT3=CT6,
  Lock3=Lock6, Reset3=Reset6, Request3=Request6,
  c31=c62, t3=t6, WD3=WD6, R3=R6, cb3=cb6,
  Trp_prm1=Trp_prm2, Lock_r1=Lock_r2, False_trip1=False_trip2,
  Reset_r1=Reset_r2, Locked_r1=Locked_r2, Trp_bkp1=Trp_bkp2,
  Arm_pR1=Arm_pR2, Arm_bR1=Arm_bR2
] endmodule

module Relay_R3 = Relay_R1 [
  R1=R3, WD1=WD3, R1_locked=R3_locked, WD5=WD1, Lock1=Lock3,
  Timer_bR1=Timer_bR3, Ttrip_pR1=Ttrip_pR3, Ttrip_bR1=Ttrip_bR3,
  FC1=FC2, CT1=CT3, FC2=FC3, CT3=CT5,
  Lock3=Lock5, Reset3=Reset5, Request3=Request5,
  c31=c53, t3=t5, WD3=WD5, R3=R5, cb3=cb5,
  Trp_prm1=Trp_prm3, Lock_r1=Lock_r3, False_trip1=False_trip3,
  Reset_r1=Reset_r3, Locked_r1=Locked_r3, Trp_bkp1=Trp_bkp3,
  Arm_pR1=Arm_pR3, Arm_bR1=Arm_bR3
] endmodule

module Relay_R4 = Relay_R1 [
  R1=R4, WD1=WD4, R1_locked=R4_locked, WD5=WD6, Lock1=Lock4,
  Timer_bR1=Timer_bR4, Ttrip_pR1=Ttrip_pR4, Ttrip_bR1=Ttrip_bR4,
  FC1=FC2, CT1=CT4, FC2=FC1, CT3=CT2,
  Lock3=Lock2, Reset3=Reset2, Request3=Request2,
  c31=c24, t3=t2, WD3=WD2, R3=R2, cb3=cb2,
  Trp_prm1=Trp_prm4, Lock_r1=Lock_r4, False_trip1=False_trip4,
  Reset_r1=Reset_r4, Locked_r1=Locked_r4, Trp_bkp1=Trp_bkp4,
  Arm_pR1=Arm_pR4, Arm_bR1=Arm_bR4
] endmodule

module Relay_R5 = Relay_R1 [
  R1=R5, WD1=WD5, R1_locked=R5_locked, WD5=WD3, Lock1=Lock5,
  Timer_bR1=Timer_bR5, Ttrip_pR1=Ttrip_pR5, Ttrip_bR1=Ttrip_bR5,
  FC1=FC3, CT1=CT5, FC2=FC1, CT3=CT1,
  Lock3=Lock1, Reset3=Reset1, Request3=Request1,
  c31=c15, t3=t1, WD3=WD1, R3=R1, cb3=cb1,
  Trp_prm1=Trp_prm5, Lock_r1=Lock_r5, False_trip1=False_trip5,
  Reset_r1=Reset_r5, Locked_r1=Locked_r5, Trp_bkp1=Trp_bkp5,
  Arm_pR1=Arm_pR5, Arm_bR1=Arm_bR5
] endmodule

module Relay_R6 = Relay_R1 [
  R1=R6, WD1=WD6, R1_locked=R6_locked, WD5=WD2, Lock1=Lock6,
  Timer_bR1=Timer_bR6, Ttrip_pR1=Ttrip_pR6, Ttrip_bR1=Ttrip_bR6,
  FC1=FC3, CT1=CT6, FC2=FC2, CT3=CT4,
  Lock3=Lock4, Reset3=Reset4, Request3=Request4,
  c31=c46, t3=t4, WD3=WD4, R3=R4, cb3=cb4,
  Trp_prm1=Trp_prm6, Lock_r1=Lock_r6, False_trip1=False_trip6,
  Reset_r1=Reset_r6, Locked_r1=Locked_r6, Trp_bkp1=Trp_bkp6,
  Arm_pR1=Arm_pR6, Arm_bR1=Arm_bR6
] endmodule

// ================================================================
// CTM modules
// ================================================================
module CTM_R1
  CT1 : [0..3] init 0;
  [] CT1 = 0 & FC1 = 1 & R1 = 5 & R5 = 5 -> (CT1' = CTM_R1);
endmodule

module CTM_R2 = CTM_R1 [CT1=CT2, FC1=FC1, R1=R2, R5=R4, CTM_R1=CTM_R2] endmodule
module CTM_R3 = CTM_R1 [CT1=CT3, FC1=FC2, R1=R3, R5=R1, CTM_R1=CTM_R3] endmodule
module CTM_R4 = CTM_R1 [CT1=CT4, FC1=FC2, R1=R4, R5=R6, CTM_R1=CTM_R4] endmodule
module CTM_R5 = CTM_R1 [CT1=CT5, FC1=FC3, R1=R5, R5=R3, CTM_R1=CTM_R5] endmodule
module CTM_R6 = CTM_R1 [CT1=CT6, FC1=FC3, R1=R6, R5=R2, CTM_R1=CTM_R6] endmodule

// ================================================================
// Signal dispatchers
// ================================================================
module Signal_Disp_R1
  Lock1   : bool init false;
  Reset1  : bool init false;
  Request1: bool init false;
  Break1  : bool init false;

  // Primary line FC1
  [Sig1] R1 = 5 & Lock1 = false & FC1 = 1 & (CT1 > 1 | WD5 = 1) ->
        (Lock1' = true);
  [Sig1] R1 = 5 & Lock1 = false & FC1 = 1 & CT1 = 1 ->
        (Lock1' = true);

  [Trp_prm1] Reset1 = false ->
        (Reset1' = true);

  [Sig1] Request1 = false & R1 = 2 ->
        (Request1' = true);

  [Sig1] Request1 = false & cb1 = 2 & R1 = 1 ->
        (Request1' = true);

  [Sig1] Break1 = false & R1 = 1 ->
        (Break1' = true);

  [Sup1] sv1 = 1 & Break1 = false ->
        (Break1' = true);

  // Backup / FC2
  [Trp_bkp1] Lock1 = false & (R1 = 4 | R1 = 5) & FC2 = 1 ->
        (Lock1' = true);

  [False_trip1] Lock1 = false & R1 = 5 & FC2 = 1 ->
        (Lock1' = true);

  [] Reset1 = false & (R1 = 1 | R1 = 2) & FC2 = 1 ->
        (Reset1' = true);
endmodule

module Signal_Disp_R2 = Signal_Disp_R1 [
  Lock1=Lock2, Reset1=Reset2, Request1=Request2, Break1=Break2,
  CT1=CT2, WD5=WD4, R1=R2, cb1=cb2, FC1=FC1, FC2=FC3, WD1=WD2,
  Sig1=Sig2, Trp_prm1=Trp_prm2, Trp_bkp1=Trp_bkp2,
  Sup1=Sup2, sv1=sv2, False_trip1=False_trip2
] endmodule

module Signal_Disp_R3 = Signal_Disp_R1 [
  Lock1=Lock3, Reset1=Reset3, Request1=Request3, Break1=Break3,
  CT1=CT3, WD5=WD1, cb1=cb3, R1=R3, FC1=FC2, FC2=FC3, WD1=WD3,
  Sig1=Sig3, Trp_prm1=Trp_prm3, Trp_bkp1=Trp_bkp3,
  Sup1=Sup3, sv1=sv3, False_trip1=False_trip3
] endmodule

module Signal_Disp_R4 = Signal_Disp_R1 [
  Lock1=Lock4, Reset1=Reset4, Request1=Request4, Break1=Break4,
  CT1=CT4, WD5=WD6, cb1=cb4, R1=R4, FC1=FC2, FC2=FC1, WD1=WD4,
  Sig1=Sig4, Trp_prm1=Trp_prm4, Trp_bkp1=Trp_bkp4,
  Sup1=Sup4, sv1=sv4, False_trip1=False_trip4
] endmodule

module Signal_Disp_R5 = Signal_Disp_R1 [
  Lock1=Lock5, Reset1=Reset5, Request1=Request5, Break1=Break5,
  CT1=CT5, WD5=WD3, cb1=cb5, R1=R5, FC1=FC3, FC2=FC1, WD1=WD5,
  Sig1=Sig5, Trp_prm1=Trp_prm5, Trp_bkp1=Trp_bkp5,
  Sup1=Sup5, sv1=sv5, False_trip1=False_trip5
] endmodule

module Signal_Disp_R6 = Signal_Disp_R1 [
  Lock1=Lock6, Reset1=Reset6, Request1=Request6, Break1=Break6,
  CT1=CT6, WD5=WD2, cb1=cb6, R1=R6, FC1=FC3, FC2=FC2, WD1=WD6,
  Sig1=Sig6, Trp_prm1=Trp_prm6, Trp_bkp1=Trp_bkp6,
  Sup1=Sup6, sv1=sv6, False_trip1=False_trip6
] endmodule

// ================================================================
// Point-to-point relay channels (Lock/Reset/Request)
// ================================================================
module Channel_R1R5
  c15 : [0..2] init 0;
  t1  : [0..4] init 0;

  [] c15 = 0 & t1 = 0 & Lock1 = true ->
       1 - COM : (c15' = 1) & (t1' = 1)
     +     COM : (c15' = 2) & (t1' = 1);

  [] c15 = 1 & R1 = 1 & R5 = 3 & Reset1 = true & t1 = 1 ->
       1 - COM : (c15' = 1) & (t1' = t1 + 1)
     +     COM : (c15' = 2) & (t1' = t1 + 1);

  [] c15 = 0 & R1 = 2 & WD1 = 1 & Request1 = true & t1 = 0 ->
       1 - COM : (c15' = 1) & (t1' = 4)
     +     COM : (c15' = 2) & (t1' = 4);

  [] c15 = 1 & R1 = 2 & WD1 = 2 & R5 = 3 &
     (Reset1 = true | Request1 = true) & t1 = 1 ->
       1 - COM : (c15' = 1) & (t1' = t1 + 1)
     +     COM : (c15' = 2) & (t1' = t1 + 1);

  [] c15 = 1 & R1 = 1 & R5 = 4 & cb1 = 2 & Request1 = true & t1 = 2 ->
       1 - COM : (c15' = 1) & (t1' = t1 + 1)
     +     COM : (c15' = 2) & (t1' = t1 + 1);
endmodule

module Channel_R2R4 = Channel_R1R5 [
  c15=c24, Request1=Request2, Lock1=Lock2,
  Reset1=Reset2, t1=t2, cb1=cb2, R1=R2, R5=R4, WD1=WD2
] endmodule

module Channel_R3R1 = Channel_R1R5 [
  c15=c31, Request1=Request3, Lock1=Lock3,
  Reset1=Reset3, t1=t3, cb1=cb3, R1=R3, R5=R1, WD1=WD3
] endmodule

module Channel_R4R6 = Channel_R1R5 [
  c15=c46, Request1=Request4, Lock1=Lock4,
  Reset1=Reset4, t1=t4, cb1=cb4, R1=R4, R5=R6, WD1=WD4
] endmodule

module Channel_R5R3 = Channel_R1R5 [
  c15=c53, Request1=Request5, Lock1=Lock5,
  Reset1=Reset5, t1=t5, cb1=cb5, R1=R5, R5=R3, WD1=WD5
] endmodule

module Channel_R6R2 = Channel_R1R5 [
  c15=c62, Request1=Request6, Lock1=Lock6,
  Reset1=Reset6, t1=t6, cb1=cb6, R1=R6, R5=R2, WD1=WD6
] endmodule

// ================================================================
// Breaker communication channels
// ================================================================
module comm_B1
  cb1 : [0..2] init 0;

  [Com_b1] cb1 = 0 & Break1 = true ->
       1 - COM : (cb1' = 1)
     +     COM : (cb1' = 2);
endmodule

module comm_B2 = comm_B1 [cb1=cb2, Break1=Break2, Com_b1=Com_b2] endmodule
module comm_B3 = comm_B1 [cb1=cb3, Break1=Break3, Com_b1=Com_b3] endmodule
module comm_B4 = comm_B1 [cb1=cb4, Break1=Break4, Com_b1=Com_b4] endmodule
module comm_B5 = comm_B1 [cb1=cb5, Break1=Break5, Com_b1=Com_b5] endmodule
module comm_B6 = comm_B1 [cb1=cb6, Break1=Break6, Com_b1=Com_b6] endmodule

// ================================================================
// Supervisory service module
// ================================================================
module Sup_serv1
  sv1 : [0..1] init 0;
  sv2 : [0..1] init 0;
  sv3 : [0..1] init 0;
  sv4 : [0..1] init 0;
  sv5 : [0..1] init 0;
  sv6 : [0..1] init 0;

  // Supervisory based on logical/timing conditions
  [a] sv1 = 0 & SUP_DUE1 -> (sv1' = 1);
  [a] sv2 = 0 & SUP_DUE2 -> (sv2' = 1);
  [a] sv3 = 0 & SUP_DUE3 -> (sv3' = 1);
  [a] sv4 = 0 & SUP_DUE4 -> (sv4' = 1);
  [a] sv5 = 0 & SUP_DUE5 -> (sv5' = 1);
  [a] sv6 = 0 & SUP_DUE6 -> (sv6' = 1);

  // Thermal hard triggers per line
  [b] sv4 = 0 & FC1 = 1 & tglobl >= T_WAIT & !SUP_DUE4 & BAD_COORD_L1_sv4 ->
      (sv4' = 1);
  [b] sv5 = 0 & FC1 = 1 & tglobl >= T_WAIT & !SUP_DUE5 & BAD_COORD_L1_sv5 ->
      (sv5' = 1);

  [b] sv1 = 0 & FC2 = 1 & tglobl >= T_WAIT & !SUP_DUE1 & BAD_COORD_L2_sv1 ->
      (sv1' = 1);
  [b] sv6 = 0 & FC2 = 1 & tglobl >= T_WAIT & !SUP_DUE6 & BAD_COORD_L2_sv6 ->
      (sv6' = 1);

  [b] sv2 = 0 & FC3 = 1 & tglobl >= T_WAIT & !SUP_DUE2 & BAD_COORD_L3_sv2 ->
      (sv2' = 1);
  [b] sv3 = 0 & FC3 = 1 & tglobl >= T_WAIT & !SUP_DUE3 & BAD_COORD_L3_sv3 ->
      (sv3' = 1);
endmodule

// ================================================================
// Rewards and labels for property specification
// ================================================================
rewards "time_to_clear"
  Fault : 1;
endrewards

label "ANY_ISOL"          = L1_Isol | L2_Isol | L3_Isol;
label "FaultPastThermal"  = Fault & (tglobl >= T_WAIT);

label "L1_Risk" = FC1 = 1 & ((sv5 = 1 & cb5 = 2) | (sv4 = 1 & cb4 = 2));
label "L2_Risk" = FC2 = 1 & ((sv1 = 1 & cb1 = 2) | (sv6 = 1 & cb6 = 2));
label "L3_Risk" = FC3 = 1 & ((sv3 = 1 & cb3 = 2) | (sv2 = 1 & cb2 = 2));

label "L1_False_Trip" =
  FC1 = 1 &
  ((CT1 = 2 & Lock1 = true & c15 = 2 & t1 = 1 & R1 = 1 & R5 = 1) |
   (CT2 = 2 & Lock2 = true & c24 = 2 & t2 = 1 & R2 = 1 & R4 = 1));

label "L2_False_Trip" =
  FC2 = 1 &
  ((CT3 = 2 & Lock3 = true & c31 = 2 & t3 = 1 & R3 = 1 & R1 = 1) |
   (CT4 = 2 & Lock4 = true & c46 = 2 & t4 = 1 & R4 = 1 & R6 = 1));

label "L3_False_Trip" =
  FC3 = 1 &
  ((CT5 = 2 & Lock5 = true & c53 = 2 & t5 = 1 & R5 = 1 & R3 = 1) |
   (CT6 = 2 & Lock6 = true & c62 = 2 & t6 = 1 & R6 = 1 & R2 = 1));

label "L1_Isolated" =
  FC1 = 1 &
  (((Break1 = true & cb1 = 1) | (Break5 = true & cb5 = 1)) &
   ((Break2 = true & cb2 = 1) | (Break4 = true & cb4 = 1)) &
   (!(sv5 = 0 & cb5 = 2) & !(sv4 = 0 & cb4 = 2)));

label "L2_Isolated" =
  FC2 = 1 &
  (((Break3 = true & cb3 = 1) | (Break1 = true & cb1 = 1)) &
   ((Break4 = true & cb4 = 1) | (Break6 = true & cb6 = 1)) &
   (!(sv1 = 0 & cb1 = 2) & !(sv6 = 0 & cb6 = 2)));

label "L3_Isolated" =
  FC3 = 1 &
  (((Break5 = true & cb5 = 1) | (Break3 = true & cb3 = 1)) &
   ((Break6 = true & cb6 = 1) | (Break2 = true & cb2 = 1)) &
   (!(sv3 = 0 & cb3 = 2) & !(sv2 = 0 & cb2 = 2)));

label "L1_Isol_Failed" = L1_Fail;

label "L2_Isol_Failed" = L2_Fail;

label "L3_Isol_Failed" = L3_Fail;
