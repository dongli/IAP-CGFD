load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

    scheme_index = 0

    available_schemes = (/ \
        (/"semi_lagrangian",         "Semi-Lagrangian scheme"/) \
    /)

    file_prefix = available_schemes(0)
    title = available_schemes(1)
    root = "."
    start_time = 0
    end_time = 125
    time_step = 1
    time_stamp = "%3.3d"

    fs = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '"+root+"/"+file_prefix+"."+time_stamp+".nc\n' $i; done")

    wks = gsn_open_wks("pdf", file_prefix)

    cmap = (/(/1.0, 1.0, 1.0 /),(/0.0, 0.0, 0.0 /), \
             (/1.00,0.91,0.91/),(/0.84,0.84,1.00/), \
             (/0.38,0.4, 0.98/),(/0.04,0.15,0.89/), \
             (/0.06,0.4, 0.46/),(/0.15,0.71,0.14/), \
             (/0.54,0.85,0.17/),(/0.98,0.98,0.22/), \
             (/0.99,0.53,0.14/),(/0.99,0.09,0.11/)/)
    gsn_define_colormap(wks, cmap)

    res = True
    res@cnLinesOn = False
    res@cnFillOn = True
    res@gsnSpreadColors = True
    res@cnLevelSelectionMode = "ManualLevels"
    res@cnMinLevelValF = -0
    res@cnMaxLevelValF = 1.0
    res@cnLevelSpacingF = 0.1
    res@lbLabelAngleF = 45
    res@tiMainString = title

    do k = 0, dimsizes(fs)-1
        f = addfile(fs(k), "r")
        system("echo 'File "+fs(k)+"' with total mass "+sprintf("%20.10f", sum(f->f)))
        plot = gsn_csm_contour(wks, f->f(0,:,:), res)
    end do

end
