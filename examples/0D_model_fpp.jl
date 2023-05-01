
using FUSE
using BoundaryPlasmaModels

ini, act = FUSE.case_parameters(:FPP; version=:v1_demount, init_from=:scalars);
dd = IMAS.dd()#
FUSE.init(dd, ini, act; do_plot=true);
actor = FUSE.ActorDivertorHeatFlux(dd, act);
BoundaryPlasmaModels.summary(actor.model)
#FUSE.ActorCXbuild(dd, act; do_plot=true, rebuild_wall=true);

#CATDEMO
dd, ini, act = FUSE.init(:CAT);
act.ActorDivertorHeatFlux.setup.lengyel.sol.f_imp = [0.02]
act.ActorDivertorHeatFlux.setup.lengyel.sol.imp = [:Ne]
actor = FUSE.ActorDivertorHeatFlux(dd, act);

BoundaryPlasmaModels.show_summary(actor.model);