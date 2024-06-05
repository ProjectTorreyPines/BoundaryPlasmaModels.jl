
using FUSE
ini, act = FUSE.case_parameters(:FPP);
dd = IMAS.dd()#
FUSE.init(dd, ini, act; do_plot=false);
#switch model
act.ActorDivertors.heat_flux_model.model = :stangeby
# check that parameters have switched Toroidal
act.ActorDivertors.heat_flux_model.setup isa FUSE.BoundaryPlasmaModels.StangebyModelParameters

actor = FUSE.ActorDivertors(dd, act);
summary(actor.boundary_plasma_models[1])
#FUSE.ActorCXbuild(dd, act; do_plot=true, rebuild_wall=true);

#CATDEMO
dd, ini, act = FUSE.init(:CAT);
act.ActorDivertors.setup.lengyel.sol.f_imp = [0.02]
act.ActorDivertorHeatFlux.setup.lengyel.sol.imp = [:Ne]
actor = FUSE.ActorDivertorHeatFlux(dd, act);

BoundaryPlasmaModels.show_summary(actor.model);