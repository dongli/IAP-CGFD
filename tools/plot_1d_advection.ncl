load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

	root = "."
	start_time = 0
	end_time = 100
	time_step = 1
	file_prefix = "upwind"
	time_stamp = "%3.3d"

	fs = systemfunc("for (( i = "+start_time+"; i <= "+end_time+"; "+ \
        "i = i+"+time_step+" )); "+ \
        "do printf '"+root+"/"+file_prefix+"."+time_stamp+".nc\n' $i; done")

	res = True
	res@vpWidthF = 0.6
	res@vpHeightF = 0.4
	res@trXMinF = 0
	res@trXMaxF = 1
	res@trYMinF = 0
	res@trYMaxF = 1
	res@tmXBTickStartF = 0
	res@tmXBTickEndF = 1
	res@tmYLTickStartF = 0
	res@tmYLTickEndF = 1
	res@xyLineThicknessF = 5.0

	wks = gsn_open_wks("pdf", file_prefix)

	do k = 0, dimsizes(fs)-1
		f = addfile(fs(k), "r")
		system("echo 'File "+fs(k)+"' with total mass "+sprintf("%20.10f", sum(f->f)))
		plot = gsn_csm_xy(wks, f->x, f->f, res)
	end do

end