using FUSE
using PlasmaFacingSurfaces
using BoundaryPlasmaModels
import BoundaryPlasmaModels
using Plots
using Logging
using LoggingExtras
#ENV["JULIA_DEBUG"] = PlasmaFacingSurfaces
#pgfplotsx()
#fig_path = "/Users/jeromeguterl/Dropbox/0D_model/design/figures/"

plotly()
pfs_logger  = TransformerLogger(global_logger()) do log
    # HTTP.jl utilizes internal modules so call parentmodule(...)
    if (log._module === PlasmaFacingSurfaces) && log.level === Logging.Debug
        # Merge can be used to construct a new NamedTuple
        # which effectively is the overwriting of fields of a NamedTuple
        return merge(log, (; level=Logging.Info))
    else
        return log
    end
end

global_logger(pfs_logger)

 ini, act = FUSE.case_parameters(:FPP; version=:v1_demount, init_from=:scalars);
 ini, act = FUSE.case_parameters(:ITER; init_from=:scalars);
 dd = IMAS.dd()#
 FUSE.init(dd, ini, act; do_plot=true);
# FUSE.IMAS.imas2json(dd,"FPP.json")
# FUSE.act2json(act,"FPP.json")
# #dd = IMAS.json2imas("PlasmaFacingSurfaces/examples/highbetap_fpp_STEP_case_with_FUSE_xpoints.json");
#dd = IMAS.json2imas("FPP.json");
#act = FUSE.json2act("FPP.json");
act.ActorPlasmaFacingSurfaces.divertors.upper.inner.l_leg = 1.0
act.ActorPlasmaFacingSurfaces.divertors.upper.outer.l_leg = 1.0
act.ActorPlasmaFacingSurfaces.divertors.lower.outer.l_leg = 1.0
act.ActorPlasmaFacingSurfaces.divertors.lower.inner.l_leg = 1.0
eqt = dd.equilibrium.time_slice[]

plt = PlasmaFacingSurfaces.plot_eq(dd, aspect_ratio=:equal, legend=:right)
xlims!(plt,PlasmaFacingSurfaces.get_xlim(plt)...)
pfs = PlasmaFacingSurfaces.PFSDesign(dd.equilibrium.time_slice[])
#PlasmaFacingSurfaces.plot_design_points()
pfs( act.ActorPlasmaFacingSurfaces; do_plot=true,do_plot_design_points=true)
plot!(pfs.divertors.lower.main_xpoint; size=(1000,1000))
x_omp = PlasmaFacingSurfaces.PlasmaFacingSurfaces.get_omp_limiter(eqt,0.002)
# PlasmaFacingSurfaces.get_l_parallel(x_omp,pfs.divertors.upper.outer.target,eqt)
# PlasmaFacingSurfaces.get_connected_mfs(x_omp,pfs.divertors.upper.outer,eqt)
get_projected_area(x_omp,λ,upper.outer)
get_connection_length(x_omp,divertor.inner)

BoundaryPlasmaModels.λq_eich(eqt)
BoundaryPlasmaModels.summary(5e6,pfs)
#FUSE.ActorCXbuild(dd, act; do_plot=true, rebuild_wall=true);