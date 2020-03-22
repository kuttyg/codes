begin

;*********** reading files************

 if(n1 .lt. 10) then
  f1="analysis.e00"+n1
  f2="analysis.e0"+n2
 else
  f1="analysis.e0"+n1
  f2="analysis.e0"+n2
 end if 
  a = addfile(f1,"r")
  b = addfile(f2,"r")

;*********** reading variables to replace **********

  lat          = wrf_user_getvar(a,"XLAT",0)
  lon          = wrf_user_getvar(a,"XLONG",0)
  s_lat_u      = wrf_user_getvar(a,"XLAT_U",0)
  s_lat_v      = wrf_user_getvar(a,"XLAT_V",0)
  s_long_u     = wrf_user_getvar(a,"XLONG_U",0)
  s_long_v     = wrf_user_getvar(a,"XLONG_V",0)
  slp2d        = wrf_user_getvar(a,"slp",0)
  dims         = dimsizes(slp2d)
  slp1d        = ndtooned(slp2d)
  smin         = minind(slp1d)
  minij        = ind_resolve(ind(slp1d.eq.min(slp2d)),dims)
  PI           = get_pi("float" )
  R            = 6371e3
  U            = a->U
  Ub           = b->U
  V            = a->V
  Vb           = b->V
  W            = a->W
  Wb           = b->W
  T            = a->T
  Tb           = b->T
  PH           = a->PH
  PHb          = b->PH
  MU           = a->MU
  MUb          = b->MU
  QVAPOR       = a->QVAPOR
  QVAPORb      = b->QVAPOR  
  
  do k1=0,348
    do k2=0,349
       phi1       = s_lat_u(k1,k2)*(PI/180)
       phi2       = s_lat_u(minij(0,0),minij(0,1))*(PI/180)
       dphi       = ((s_lat_u(k1,k2)-s_lat_u(minij(0,0),minij(0,1))))*(PI/180)
       dsci       = ((s_long_u(k1,k2)-s_long_u(minij(0,0),minij(0,1))))*(PI/180)
       x          = (sin(dphi/2))^2+cos(phi1)*cos(phi2)*(sin(dsci/2))^2
       c          = 2*atan2(sqrt(x),sqrt(1-x))
       dis        = (R*c)/1000
       if(dis.lt.300) then
          alpha = 1.0
          U(:,:,k1,k2) = alpha*U(:,:,k1,k2) + (1-alpha)* Ub(:,:,k1,k2)           
       else if(dis.lt.600 .and. dis.gt.300) then
          d     = (dis-300)
          alpha = .42 + 0.5*cos((PI*d)/300) + 0.08*cos((2*PI*d)/300)
          U(:,:,k1,k2) = alpha*U(:,:,k1,k2) + (1-alpha)* Ub(:,:,k1,k2)           
       else
          alpha = 0.0
          U(:,:,k1,k2) = alpha*U(:,:,k1,k2) + (1-alpha)* Ub(:,:,k1,k2)
       end if
       end if
    end do
  end do
  do k1=0,349
    do k2=0,348
       phi1       = s_lat_v(k1,k2)*(PI/180)
       phi2       = s_lat_v(minij(0,0),minij(0,1))*(PI/180)
       dphi       = ((s_lat_v(k1,k2)-s_lat_v(minij(0,0),minij(0,1))))*(PI/180)
       dsci       = ((s_long_v(k1,k2)-s_long_v(minij(0,0),minij(0,1))))*(PI/180)
       x          = (sin(dphi/2))^2+cos(phi1)*cos(phi2)*(sin(dsci/2))^2
       c          = 2*atan2(sqrt(x),sqrt(1-x))
       dis        = (R*c)/1000
       if(dis.lt.300) then
          alpha = 1.0
          V(:,:,k1,k2) = alpha*V(:,:,k1,k2) + (1-alpha)* Vb(:,:,k1,k2)           
       else if(dis.lt.600 .and. dis.gt.300) then
          d = (dis-300)
          alpha = .42 + 0.5*cos((PI*d)/300) + 0.08*cos((2*PI*d)/300)
          V(:,:,k1,k2) = alpha*V(:,:,k1,k2) + (1-alpha)* Vb(:,:,k1,k2)           
       else
          alpha = 0.0
          V(:,:,k1,k2) = alpha*V(:,:,k1,k2) + (1-alpha)* Vb(:,:,k1,k2)
       end if
       end if
    end do
  end do
  do k1=0,348
    do k2=0,348
       phi1       = lat(k1,k2)*(PI/180)
       phi2       = lat(minij(0,0),minij(0,1))*(PI/180)
       dphi       = ((lat(k1,k2)-lat(minij(0,0),minij(0,1))))*(PI/180)
       dsci       = ((lon(k1,k2)-lon(minij(0,0),minij(0,1))))*(PI/180)
       x          = (sin(dphi/2))^2+cos(phi1)*cos(phi2)*(sin(dsci/2))^2
       c          = 2*atan2(sqrt(x),sqrt(1-x))
       dis        = (R*c)/1000
       if(dis.lt.300) then
          alpha = 1.0
          W(:,:,k1,k2) = alpha*W(:,:,k1,k2) + (1-alpha)* Wb(:,:,k1,k2)
          T(:,:,k1,k2) = alpha*T(:,:,k1,k2) + (1-alpha)* Tb(:,:,k1,k2)
          PH(:,:,k1,k2) = alpha*PH(:,:,k1,k2) + (1-alpha)* PHb(:,:,k1,k2)
          QVAPOR(:,:,k1,k2) = alpha*QVAPOR(:,:,k1,k2) + (1-alpha)* QVAPORb(:,:,k1,k2)
          MU(:,k1,k2) = alpha*MU(:,k1,k2) + (1-alpha)* MUb(:,k1,k2)                     
       else if(dis.lt.600 .and. dis.gt.300) then
          d = (dis-300)
          alpha = .42 + 0.5*cos((PI*d)/300) + 0.08*cos((2*PI*d)/300)
          W(:,:,k1,k2) = alpha*W(:,:,k1,k2) + (1-alpha)* Wb(:,:,k1,k2)
          T(:,:,k1,k2) = alpha*T(:,:,k1,k2) + (1-alpha)* Tb(:,:,k1,k2)
          PH(:,:,k1,k2) = alpha*PH(:,:,k1,k2) + (1-alpha)* PHb(:,:,k1,k2)
          QVAPOR(:,:,k1,k2) = alpha*QVAPOR(:,:,k1,k2) + (1-alpha)* QVAPORb(:,:,k1,k2)
          MU(:,k1,k2) = alpha*MU(:,k1,k2) + (1-alpha)* MUb(:,k1,k2)           
       else
          alpha = 0.0
          W(:,:,k1,k2) = alpha*W(:,:,k1,k2) + (1-alpha)* Wb(:,:,k1,k2)
          T(:,:,k1,k2) = alpha*T(:,:,k1,k2) + (1-alpha)* Tb(:,:,k1,k2)
          PH(:,:,k1,k2) = alpha*PH(:,:,k1,k2) + (1-alpha)* PHb(:,:,k1,k2)
          QVAPOR(:,:,k1,k2) = alpha*QVAPOR(:,:,k1,k2) + (1-alpha)* QVAPORb(:,:,k1,k2)
          MU(:,k1,k2) = alpha*MU(:,k1,k2) + (1-alpha)* MUb(:,k1,k2)
       end if
       end if
    end do
  end do
  f = addfile(f1,"w")
  f->U = U
  f->V = V
  f->W = W
  f->T = T
  f->QVAPOR = QVAPOR
  f->PH = PH
  f->MU = MU
end
