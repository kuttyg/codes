begin

;*********** reading files************
  a = addfile("analysis.e012","r")

;*********** reading variables to replace **********

  lat          = wrf_user_getvar(a,"XLAT",0)
  lon          = wrf_user_getvar(a,"XLONG",0)
  T            = a->T  
  A            = 1.5
; ********Applying the formula on all relevant grids*********
  do k1=0,348
    do k2=0,348
      if (lat(k1,k2) .gt. 10)
          delT = (lat(k1,k2)-10)/10
          T(:,:,k1,k2)=T(:,:,k1,k2)+A*delT
      end if
    end do
  end do
  f = addfile("analysis.e012","w")
  f->T = T
end
