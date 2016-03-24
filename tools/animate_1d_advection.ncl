load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    available_schemes = (/ \
        (/"upwind",         "Upwind scheme"/), \
        (/"beam_warming",   "Beam-Warming scheme"/), \
        (/"lax_wendroff",   "Lax-Wendroff scheme"/), \
        (/"fromm",          "Fromm scheme"/), \
        (/"ftcs",           "Forward-time center-space scheme"/), \
        (/"leapfrog",       "Leap-frog scheme"/), \
        (/"crank_nicolson", "Crank-Nicolson scheme"/), \
        (/"mpdata",         "MPDATA scheme"/), \
        (/"tspas",          "TSPAS scheme"/), \
        (/"ffsl",           "FFSL scheme"/) \
    /)

    if (isdefined("scheme")) then
        do i = 0, 9
            if (available_schemes(i,0) .eq. scheme) then
                scheme_index = i
                break
            end if
        end do
    else
        scheme_index = 8
    end if

    file_prefix = available_schemes(scheme_index, 0)
    title = available_schemes(scheme_index, 1)
    root = "."
    start_time = 0
    if (.not. isdefined("end_time")) then
        end_time = 126
    end if
    if (.not. isdefined("time_step")) then
        time_step = 1
    end if
    time_stamp = "%3.3d"

    fs = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '"+root+"/"+file_prefix+"."+time_stamp+".nc\n' $i; done")

    res = True
    res@vpWidthF = 0.6
    res@vpHeightF = 0.4
    res@trXMinF = 0
    res@trXMaxF = 1
    res@trYMinF = -0.2
    res@trYMaxF = 1.2
    res@tmXBTickStartF = res@trXMaxF
    res@tmXBTickEndF = res@trXMinF
    res@tmYLTickStartF = res@trYMinF
    res@tmYLTickEndF = res@trYMaxF
    res@xyLineThicknessF = 5.0
    res@tiMainString = title

    wks = gsn_open_wks("pdf", file_prefix)

    do k = 0, dimsizes(fs)-1
        f = addfile(fs(k), "r")
        system("echo 'File "+fs(k)+"' with total mass "+sprintf("%20.10f", sum(f->f)))
        plot = gsn_csm_xy(wks, f->x, f->f, res)
    end do

end
