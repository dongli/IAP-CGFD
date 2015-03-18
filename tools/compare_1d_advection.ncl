load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin

	file_paths = (/"mpdata.200.iord1.nc", \
	               "mpdata.200.iord2.nc", \
	               "mpdata.200.iord3.nc", \
	               "mpdata.200.iord4.nc"/)

	colors = (/"black","red","cyan","blue"/)

	wks = gsn_open_wks("pdf", "compare")

    res = True
    res@gsnDraw = False
    res@gsnFrame = False
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

    plots = new(dimsizes(file_paths), graphic)

    do i = 0, dimsizes(plots)-1
    	f = addfile(file_paths(i), "r")
    	res@xyLineColor = colors(i)
    	plots(i) = gsn_csm_xy(wks, f->x, f->f, res)
    end do

    do i = 1, dimsizes(plots)-1
	    overlay(plots(0), plots(i))
	end do
	draw(plots(0))
	frame(wks)

end