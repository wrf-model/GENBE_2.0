#! /bin/ksh
#-----------------------------------------------------------------------
# Script gen_be_wrapper.ksh
#
# Purpose: Calculates background error statistics for WRF-Var.
#-----------------------------------------------------------------------

#[1] Define job by overriding default environment variables:
# this steps need to be run successively by setting up each variable to true

export RUN_GEN_BE_STAGE0=true  # Run stage 0 (create perturbation files).
export RUN_GEN_BE_STAGE1=true   # Run stage 1 (Remove mean, split variables).
export RUN_GEN_BE_STAGE2=true  # Run stage 2 (Regression coefficients).
export RUN_GEN_BE_STAGE3=true  # Run stage 3 (Vertical Covariances).
export RUN_GEN_BE_STAGE4=true  # Run stage 4 (Horizontal Covariances).
export RUN_GEN_BE_DIAGS=true   # Generate the be.nc file


export GEN_BE_DIR="/glade/p/mmm/liuz/ys_code/GENBE_2.0"  # code directory
export START_DATE=2016071900                           # the first perturbation valid date
export END_DATE=2016072100                             # the last perturbation valid date
export NUM_LEVELS=49                                   # number levels - 1
export FC_DIR="/glade/p/mmm/liuz/ys_code/exp-new_exp-ahi-1"  # forecast directory"
export RUN_DIR="/glade/p/mmm/liuz/ys_code/GENBE_2.0/scripts"            # scripts directory
export WORK_DIR="/glade/p/mmm/liuz/ys_code/exp-new_exp-ahi-1/be_lsmethod2" # working directory
export DOMAIN=01
export INTERVAL=12
export STRIDE=1
export FCST_RANGE=24
export FCST_RANGE1=24
export FCST_RANGE2=12

export BE_METHOD="NMC"  #"ENS"
export NE="1"

#[2] Run gen_be:
cd "$GEN_BE_DIR/scripts"
./gen_be.ksh
