set (wavcav_SOURCES 
   PL_cubic.c
   PL_quad_vil.c
   PL_view_vil.c
   arcplane.c
   baryc.c
   baryc_vil.c
   bezcoons.c
   bivariate.c
   circle.c
   circmesh.c
   circular.c
   circular_arc.c
   collapse.c
   composite_curve.c
   concave1.c
   concave2.c
   concave3.c
   conic.c
   consistency.c
   convertion.c
   convex1.c
   convex2.c
   convhull.c
   copy_ent.c
   deltria.c
   dijkspat.c
   drawinter.c
   drawres.c
   drawsph.c
   earscut.c
   elem_vil.c
   exinterp.c
   exinterp2.c
   expand.c
   fairing.c
   flip_vil.c
   flowcurve.c
   fuse_vil.c
   gap.c
   graphmst.c
   holefill.c
   insert.c
   intermediate.c
   line_entity.c
   loadquad_vil.c
   lod_vil.c
   mapmanif.c
   memory.c
   merge_vil.c
   misc.c
   misc_vil.c
   molec.c
   mst_vil.c
   norm_nurbs.c
   nurbs.c
   opersingle.c
   opfloater.c
   partcas1.c
   partcas2.c
   partcas3.c
   particulier.c
   planeeq_vil.c
   pln_sph.c
   polygonal.c
   preremov_vil.c
   procinterpol.c
   prw_arr.c
   quad3D_vil.c
   quadric_vil.c
   quadtri_vil.c
   qualquad.c
   rectif.c
   refine1.c
   refine2.c
   reorg_vil.c
   segm_segm.c
   selfinters1.c
   selfinters2.c
   shift_vil.c
   skite_vil.c
   smooth.c
   smoothingpln.c
   solvent.c
   solvers.c
   solvers_vil.c
   sphere.c
   submanif.c
   test.c
   test_vil.c
   toroidal1.c
   toroidal2.c
   toroidal3.c
   toroidal4.c
   toroidal5.c
   toroidpos.c
   traverse.c
   treatblend1.c
   treatblend2.c
   treatblend3.c
   treatblend4.c
   treatblend5.c
   treatblend6.c
   treatblend7.c
   treatblend8.c
   treatmanif.c
   triangulate1.c
   triangulate2.c
   trimmed_sphere.c
   trimmed_surface.c
   trimsphere1.c
   trimsphere2.c
   upd_discr1.c
   upd_discr2.c
   verticesview_vil.c
   weights_vil.c
  )

add_library(wavcav ${wavcav_SOURCES})

install(FILES pln_sph.h DESTINATION include)
install(FILES geodesic.h DESTINATION include)
install(FILES coarsequad.h DESTINATION include)
install(FILES cavity.h DESTINATION include)
install(FILES triang.h DESTINATION include)
install(FILES splinemol.h DESTINATION include)
install(FILES smooth.h DESTINATION include)
install(FILES partcas.h DESTINATION include)
install(FILES meshsas.h DESTINATION include)
install(FILES eval.h DESTINATION include)
install(FILES sas.h DESTINATION include)

install(TARGETS wavcav
	ARCHIVE DESTINATION lib)

