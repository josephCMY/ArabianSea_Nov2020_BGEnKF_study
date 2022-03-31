












program da_wrfvar_main

   
   
   
   
   
   
   
   

   use module_symbols_util, only : wrfu_finalize

   use da_control, only : trace_use, var4d
   use da_tracing, only : da_trace_init, da_trace_report, da_trace_entry, &
      da_trace_exit
   use da_wrf_interfaces, only : wrf_shutdown, wrf_message, disable_quilting
   use da_wrfvar_top, only : da_wrfvar_init1,da_wrfvar_init2,da_wrfvar_run, &
      da_wrfvar_finalize

   implicit none

   

   call disable_quilting

   call da_wrfvar_init1

   if (trace_use) call da_trace_init
   if (trace_use) call da_trace_entry("da_wrfvar_main")

   call da_wrfvar_init2

   call da_wrfvar_run

   call da_wrfvar_finalize


   call wrf_message("*** WRF-Var completed successfully ***")

   if (trace_use) call da_trace_exit("da_wrfvar_main")
   if (trace_use) call da_trace_report

   call wrfu_finalize
   call wrf_shutdown

end program da_wrfvar_main

