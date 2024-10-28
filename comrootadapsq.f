        implicit real *8 (a-h,o-z)
        dimension w(10 000 000),x(20 00),errors(10 000)
        complex *16 c0,ima
        complex *16 roots(10 000),centers(10 000 000)
        save
        data ima/(0.0d0,1.0d0)/

        external test_fun1
        external test_fun2
        external test_fun3
        external test_fun4
        external test_fun5
        external test_fun6
        external test_fun7
        external test_fun8
        external test_fun9

        call prini(6,13)


        print *,'ENTER n'
        read *,n


ccc     order of expansion
        norder=30
        norder=20
        norder=80
        norder=40

        print *,norder



        c0=0
        sqw=2.13

        c0=10-20*ima
        sqw=101

        c0=-0.02+0.01d0*ima
        c0=0
        sqw=2.13d0
        sqw=2.5



        print *,'finding root in a square with center'
        print *, c0
        print *,'and side length'
        print *, sqw

c
c if derivative values are available
c
        ifdval=1
c 
c if to use residue to count roots
c
        ifres=1
c 
c if newton refinements
c
        ifnewton=1

c 
c if extradepth to go after conv. acc.
c
        nexdepth=1


        epsexp=1d-14
        epsexp=1d-2
        epsexp=1d-4
        call mach_zero(epsexp)
        epsexp=1d-30
        epsexp=1d-12
        epsexp=1d-10

        call complex_root_adap_square(test_fun7,norder,epsexp,
     1          roots,errors,nrtot,centers,nc,c0,sqw,w,
     2          ifdval,ifres,ifnewton,nexdepth)

c
c plot all complex roots
c
 
        iw=61
        itype=2
        call zquaplot(iw,roots,nrtot,itype,
     1      'sin(50/(phi z-2)*')

        do 1000 i=1,1000
        x(i)=i
 1000 continue

        iw=62
        itype=2
        call quagraph(iw,x,log10(errors),nrtot,
     1      itype,'Log scale roots error*')

 
        iw=63
        itype=2
        call zquaplot(iw,centers,nc,itype,'All square centers*')


        call prin2('all roots found *',roots,2*nrtot)
ccc        call prin2_long('all roots found *',roots,2*nrtot)
c        call prin2('root error is *',errors,nrtot)
        call prinf('total root number is *',nrtot,1)

        stop
        end


        subroutine test_fun1(val,dval,x)
        implicit real *8 (a-h,o-z)
        complex *16 val,dval,x
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/
        
        done=1

        val=(x-0.5d0)**5*(x-0.9d0)**4*(x+0.8d0)
     1           *(x-0.7d0*ima)*(x+0.1d0*ima)**3
ccc        val=(x-0.5)*(x-0.9)**8*(x+0.8)*(x-0.7*ima)**6*(x+0.1*ima)**3

        dval=5d0*(x-0.5d0)**4*(x-0.9d0)**4*(x+0.8d0)
     1           *(x-0.7d0*ima)*(x+0.1d0*ima)**3
        dval=dval+4d0*(x-0.5d0)**5*(x+0.8d0)*(x-0.9d0)
     1         **3*(x-0.7d0*ima)*(x+0.1d0*ima)**3
        dval=dval+ (x-0.5d0)**5*(x-0.9d0)**4
     1          *(x-0.7d0*ima)*(x+0.1d0*ima)**3
        dval=dval+ (x-0.5d0)**5*(x-0.9d0)**4
     1              *(x+0.8d0)*(x+0.1d0*ima)**3
        dval=dval+3d0*(x-0.5d0)**5*(x-0.9d0)**4*(x+0.8d0)
     1      *(x-0.7d0*ima)*(x+0.1d0*ima)**2 


        return
        end


        subroutine test_fun2(val,dval,x)
        implicit real *8 (a-h,o-z)
        complex *16 val,dval,x
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/
        
        done=1
        val=(x-0.5d0)*(x-0.9d0)*(x+0.8d0)*
     1          (x-0.7d0*ima)*(x+0.1d0*ima)


        dval=(x-0.9d0)*(x+0.8d0)*(x-0.7d0*ima)*(x+0.1d0*ima)
        dval=dval+(x-0.5d0)*(x+0.8d0)*(x-0.7d0*ima)*(x+0.1d0*ima)
        dval=dval+ (x-0.5d0)*(x-0.9d0)*(x-0.7d0*ima)*(x+0.1d0*ima)
        dval=dval+ (x-0.5d0)*(x-0.9d0)*(x+0.8d0)*(x+0.1d0*ima)
        dval=dval+ (x-0.5d0)*(x-0.9d0)*(x+0.8d0)*(x-0.7d0*ima)


        return
        end



        subroutine test_fun3(val,dval,x)
        implicit real *8 (a-h,o-z)
        complex *16 val,dval,x
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/
        
        done=1
        const=1
        pi=atan(done)*4
        val=sin(x*const*3*pi/2)/(x-done*2)

        dval=const*3*pi/2*cos(x*const*3*pi/2)/(x-done*2)
        dval=dval-sin(x*const*3*pi/2)/(x-done*2)**2


        return
        end

        subroutine test_fun4(val,dval,x)
        implicit real *8 (a-h,o-z)
        complex *16 val,dval,x
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/
        
        done=1
        pi=atan(done)*4
        val=sinh(x*3*pi/2)/(x-done*2)

        dval=3*pi/2*cosh(x*3*pi/2)/(x-done*2)
        dval=dval-sinh(x*3*pi/2)/(x-done*2)**2


        return
        end


        subroutine test_fun5(val,dval,x)
        implicit real *8 (a-h,o-z)
        complex *16 val,dval,x
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/
        
        done=1
        pi=atan(done)*4
        val=sin(x*3*pi)/(x-done*2)

        dval=3*pi*cos(x*3*pi)/(x-done*2)
        dval=dval-sin(x*3*pi)/(x-done*2)**2


        return
        end

        subroutine test_fun6(val,dval,x)
        implicit real *8 (a-h,o-z)
        complex *16 val,dval,x,r0
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/
        
        done=1
        d=100
ccc        x0=1.3
        x0=2
        r0=exp(0.1*ima)
        r0=1
        pi=atan(done)*4
        val=sin(d/(x*r0-x0*done))

        dval=-cos(d/(x*r0-x0*done))*(d*r0/(x*r0-x0*done)**2)


        return
        end


        subroutine test_fun7(val,dval,x)
        implicit real *8 (a-h,o-z)
        complex *16 val,dval,x,r0
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/
        
        done=1
        d=100
ccc        x0=1.3
        x0=2
        r0=exp(0.4*ima)
ccc        r0=exp(0.1*ima*x)
        pi=atan(done)*4
        val=sin(d/(x*r0-x0*done))

        dval=-cos(d/(x*r0-x0*done))*(d*r0/(x*r0-x0*done)**2)


        return
        end

        subroutine test_fun8(val,dval,x)
        implicit real *8 (a-h,o-z)
        complex *16 val,dval,x,r0
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/
        
        done=1
        pi=atan(done)*4
        d=100
ccc        x0=1.3
        x0=2
ccc        r0=exp(0.5*ima*x)
ccc        r0=exp(0.1*ima*x-ima*pi/4)
        r0=exp(ima*pi/4)
        val=sin(d/(x*r0-x0*done))

        dval=-cos(d/(x*r0-x0*done))*(d*r0/(x*r0-x0*done)**2)


        return
        end


        subroutine test_fun9(val,dval,x)
        implicit real *8 (a-h,o-z)
        complex *16 val,dval,x,r0
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/
        
        done=1
        pi=atan(done)*4
        d=50
        phi=1.1
ccc        x0=1.3
        x0=2
ccc        r0=exp(1*ima*x)
ccc        r0=exp(1*ima*x)
ccc        r0=exp(0.5*ima*x-ima*pi/4)
        r0=exp(phi*ima*x-ima*pi/4)
        val=sin(d/(x*r0-x0*done))

        dval=-cos(d/(x*r0-x0*done))*(d/(x*r0-x0*done)**2)*
     1      (r0*ima*phi*x+r0)


        return
        end

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c        this is the end of the debugging code and the beginning
c        of the complex rootfinding subroutines
c 
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine complex_root_adap_square(fun,norder,epsexp,roots,
     1      errors,nrtot,centers,nc,c0,sqw,w,
     2      ifdval,ifres,ifnewton,nexdepth)
        implicit real *8 (a-h,o-z)
        complex *16 c0,roots(1),centers(1)
        dimension w(1)

c       This subroutine finds all complex roots of the user provided 
c       analytic function fun in the square with length sqw centered 
c       at c0. It returns all roots found in roots with nrtot the number of
c       roots. This subroutine is only a memory management subroutine;
c       all actually work is done in the subroutine complex_root_adap_square0.
c
c       The root finder proceeds by first expanding fun on the boundary 
c       of the square by preconstructed polynomials of order norder,
c       statisfies some three-term recursions.  If the expansion doesn't 
c       converge by the pre-specified exit condition, it keeps divide 
c       the square until the expansion converges. For each convergent expansion, 
c       the root finder finds all roots within the square by finding 
c       eigenvalues of the corresponding colleague matrix. After all 
c       convergent squared being solved, it removes potentially
c       duplicated roots from adjacent squares. 
c
c       The subroutine provides the option to refine all roots with
c       Newton steps. It also provides the option to estimate the 
c       number of roots inside of the square by the residue formula.
c       
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                   input parameters:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       fun: the function whose roots to be found
c
c       norder: the order of the expansion, no smaller than norder=30
c                  recommended value: norder=40 for double precision
c                                     norder=80 for extended precision
c       epsexp: the relative accuracy of polynomial expansion.
c               It sets the LOWERBOUND on the precison of 
c               the computed roots.
c
c          recommended value: from 1d-4 to 1d-12 
c                              (double precision)
c
c
c       c0: the center the square
c       sqw: length of the square
c
c       ifdval: flag for if the derivative is available
c       ifnewton: flag for newton refinements 
c       ifres: flag for residue root counting by contour int.
c       The last two flags are only meaningful if ifdval=1

c       nexdepth: go nexdepth level deeper after reaching the accuracy
c to ensure convergence
        
c
c       The the input function should be of the form
c                      fun(val,dval,z),
c       where val=f(z), the value of fun at z, and dval= df/dz (z),
c       the derivative of fun at z.
c       If the derivative value dval is not available, set dval to be
c       any value.
c
c                                WARNING:
c       The subroutine has taken precoutions to ensure the relative
c       expansion accuracy on the boundary exceeds the one set by epsexp. 
c       As a result, for eps~1d-12, it achieves machine eps precision
c       whenever it is possible,
c       so epsexp<1d-12 is in general meaningless and will only
c       prolong unnecessarily the running time. 
c
c       The parameter epsexp only sets the LOWERBOUND of the relative error 
c       of the roots: loss of 2-3 digits should be expected since the polynomial
c       approximation is on the boundary, where the function values are  
c       larger than those inside. 
c       When high accuracy is needed (close to machine eps), it is
c       recommended to turn newton refinement on (ifnewton=1).
c       Although not recommended,  using smaller
c       expansion order norder, such as 
c                              norder=15-20 for double precision 
c                              norder=30-40 for extended precision 
c       also increases the accuracy.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                   output parameters:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       roots: all roots found
c       errors: errors of all roots
c       nrtot: the total number of roots found
c              also the length of errors
c       centers: centers of all subdivided squares 
c       nc: the length of centers
c
c
c       ier: error flag
c            0, successful
c            1024, maxlevel 20 reached
c            2048, max num of squares  10 000 000 reached

c                   work array: 
c           w: at least 4330*norder+6363 long 
c


c
c memory allocation
c

c
c   For some reason, the precomputed SVD cannot be
c   loaded for norder<30, so set the minimum to be 30
c
        if(norder.lt.30) nnorder=30
        if(norder.ge.30) nnorder=norder


ccc        npt=50*4
        npt=60*4

        iz=1
        lz=npt*2+10
       
        icw=iz+lz
        lcw=npt*2+10

        iq=icw+lcw 
        lq=(norder+1)*npt*2

        iwk=iq+lq
        lwk=2*(npt*(nnorder+1)*8+500)+10

        iv=iwk+lwk
        lv=norder*2+10
        
        ivo=iv+lv
        lvo=norder*2+10

        ic=ivo+lvo
        lc=(norder+1)*2+10

        icr=ic+lc
        lcr=norder*2+10

        ir=icr+lcr
        lr=norder*2+10


ccc        call prinf('length is *',lr+ir,1)

        call complex_root_adap_square0(fun,norder,epsexp,roots,errors,
     1   nrtot,centers,nc,c0,sqw,w(iz),w(icw),w(iq),w(iwk),
     2      w(iv),w(ivo),w(ic),w(icr),w(ir),ifdval,ifres,ifnewton,
     3      nexdepth)

        return
        end
        


        subroutine complex_root_adap_square0(fun,norder,epsexp,
     1  rootslist,errlist,nrtot,c0list,ntot,c0,sqw,z,cw,q,work,vd,voffd,
     2          coefs,croots,roots,ifdval,ifres,ifnewton,nexdepth)
        implicit real *8 (a-h,o-z)
        dimension list(10,10 000 000),errlist(100 000),x(100 000),
     1              ilevel(2,25)
        complex *16 geomlist(2,10 000 000),rootslist(1),c0list(1)
        complex *16 vd(1),voffd(1),coefs(1),q(1),z(1),work(1),
     1              cw(1),croots(1),roots(1)
        complex *16 ima,c0,c0div(4),uc,uh,uv,z0,val,dval
        data ima/(0.0d0,1.0d0)/

c       This subroutine finds all complex roots of the user provided 
c       analytic function fun in the square with length sqw centered 
c       at c0. It returns all roots found in roots with nrtot the number of
c       roots. 
c
c       The root finder proceeds by first expanding fun on the boundary 
c       of the square by preconstructed polynomials of order norder,
c       statisfies some three-term recursions.  If the expansion doesn't 
c       converge by the pre-specified exit condition, it keeps divide 
c       the square until the expansion converges. For each convergent expansion, 
c       the root finder finds all roots within the square by finding 
c       eigenvalues of the corresponding colleague matrix. After all 
c       convergent squared being solved, it removes potentially
c       duplicated roots from adjacent squares. 
c
c       The subroutine provides the option to refine all roots with
c       Newton steps. It also provides the option to estimate the 
c       number of roots inside of the square by the residue formula.
c       
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                   input parameters:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c 
c       fun: the function whose roots to be found
c
c       norder: the order of the expansion
c                  recommended value: norder=40 for double precision
c                                     norder=80 for quadrupole precision
c       c0: the center the square
c       sqw: length of the square
c
c
c       epsexp: the relative accuracy of polynomial expansion.
c               It sets the LOWERBOUND on the precison of 
c               the computed roots.
c
c          recommended value: from 1d-4 to 1d-12 
c                              (double precision)

c
c       ifdval: flag for if the derivative is available
c       ifnewton: flag for newton refinements 
c       ifres: flag for residue root counting by contour int.
c       The last two flags are only meaningful if ifdval=1

c        nexdepth: go nexdepth level deeper after reaching the accuracy
c to ensure convergence
c
c       The the input function should be of the form
c                      fun(val,dval,z),
c       where val=f(z), the value of fun at z, and dval= df/dz (z),
c       the derivative of fun at z.
c       If the derivative value dval is not available, set dval to be
c       any value.
c
c                                WARNING:
c       The subroutine has taken precoutions to ensure the relative
c       expansion accuracy on the boundary exceeds the one set by epsexp. 
c       As a result, for eps~1d-12, it achieves machine eps precision
c       whenever it is possible,
c       so epsexp<1d-12 is in general meaningless and will only
c       prolong unnecessarily the running time. 
c
c       The parameter epsexp only sets the LOWERBOUND of the relative error 
c       of the roots: loss of 2-3 digits should be expected since the polynomial
c       approximation is on the boundary, where the function values are  
c       larger than those inside. 
c       When high accuracy is needed (close to machine eps), it is
c       recommended to turn newton refinement on (ifnewton=1).
c       Although not recommended,  using smaller
c       expansion order norder, such as 
c                              norder=15-20 for double precision 
c                              norder=30-40 for extended precision 
c       also increases the accuracy.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                   output parameters:
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       roots: all roots found
c       errors: errors of all roots
c       nrtot: the total number of roots found
c              also the length of errors
c       centers: centers of all subdivided squares 
c       nc: the length of centers
c       w: the array containing all precomputed arrays,
c          quad tree information
c          must be 4*50 .........long
c
c       ier: error flag
c            0, successful
c            1024, maxlevel 20 reached
c            2048, max num of squares  10 000 000 reached
c
c
c
c           work array:
c       list: information of the quadtree of divided squares
c       list(1): level, list(2): parent index, 
c       list(3):ifdone, list(4-7): children index            
c       list(10):ierqr, error ind from cqr for the box
c       geomlist: center and side length of all squares
c
c
c
c
        if(norder.lt.30) nnorder=30
        if(norder.ge.30) nnorder=norder

c
c go nexdepth level deeper after reaching the accuracy
c to ensure convergence
c

ccc        nexdepth=1
ccc        call prinf('Extra depth to go is set to be *',nexdepth,1)

c
c load orthogonal poly and three term recursion
c
        call orth_poly_load(q,vd,voffd,norder)

c
c points on the square and weights
c
        n=60
        call get_pts_and_wts(z,cw,n)
ccc        call weight_load(cw,n)

c
c get least square sol matrices
c

c
c   For some reason, the precomputed SVD cannot be
c   loaded for norder<30, so set the minimum to be 30
c
        ml=4*n*(nnorder+1)*8+500
        call work_array_load(work,ml)

c
c get fekete pts
c
ccc        mm=400
ccc        call get_fekete_pts(zfg,mm)

        print *,'loading finished'



        geomlist(1,1)=c0
        geomlist(2,1)=sqw
        
        
        level=0
        ifirst=1
        ilast=1
        ilastnew=ilast
ccc        icount=1
        nrtot=0
        ier=0
c
c level
c
        list(1,1)=level
c
c parent
c
        list(2,1)=0

        do 2000 level=0,25

        iffini=1
        call prinf('level is *',level,1)
c        call prinf('ifirst is *',ifirst,1)
c        call prinf('ilast is *',ilast,1)

        do 1000 icount=ifirst,ilast
c        call prinf('icount is *',icount,1)
        c0=geomlist(1,icount)
        sqw=real(geomlist(2,icount))

c        call prin2('c0 is *',c0,2)
c        call prin2('sqw is *',sqw,1)

        if(icount.eq.1) ifdone=2+nexdepth

        if(icount.ne.1) then
         iparent=list(2,icount)
         ifdone=list(3,iparent)
        endif
        
cccc        if(icount.eq.6) then
cccc        call prinf('icount is *',icount,1)
cccc        call prinf('iparent is *',iparent,1)
cccc        call prinf('ifdone is *',ifdone,1)
cccccc        stop
cccc        endif

        call complex_root_one_square_adap(epsexp,ifdone,fun,c0,sqw,
     1      vd,voffd,q,z,cw,n,work,croots,norder,roots,nr,ierqr)

ccc        call prinf('ifdone is *',ifdone,1)

        ifcv=0
        if(ifdone.eq.1) ifcv=1

        iffini=ifcv*iffini
        list(3,icount)=ifdone

        if(ifcv.eq.1) then
         do 1100 i=1,nr
         rootslist(i+nrtot)=roots(i)

ccc         errlist(i+nrtot)=abs(err(i))
 1100 continue
        
         nrtot=nrtot+nr


ccc         call prin2('roots are *', roots,nr*2)

        endif
c
c error ind from cqr
c

        list(10,icount)=ierqr
ccc        call prinf('icount is *',icount,1)
ccc        call prinf('ierqr at icount is *',ierqr,1)
        
        if(ifcv.eq.0) then
         call divide_square(c0,sqw,c0div,sqwdiv)

         n1=ilastnew+1
         n2=ilastnew+2
         n3=ilastnew+3
         n4=ilastnew+4
         list(1,n1)=level+1
         list(2,n1)=icount
         geomlist(1,n1)=c0div(1)
         geomlist(2,n1)=sqwdiv
         list(1,n2)=level+1
         list(2,n2)=icount
         geomlist(1,n2)=c0div(2)
         geomlist(2,n2)=sqwdiv
         list(1,n3)=level+1
         list(2,n3)=icount
         geomlist(1,n3)=c0div(3)
         geomlist(2,n3)=sqwdiv
         list(1,n4)=level+1
         list(2,n4)=icount
         geomlist(1,n4)=c0div(4)
         geomlist(2,n4)=sqwdiv

c
c children
c
         list(4,icount)=n1
         list(5,icount)=n2
         list(6,icount)=n3
         list(7,icount)=n4

         ilastnew=ilastnew+4
         
         if(ilastnew.ge.10 000 000) then
         ier=2048
         goto 3000
         endif

        endif
 1000 continue

c
c record index of each level
c

        ilevel(1,level+1)=ifirst
        ilevel(2,level+1)=ilast

        call prinf('num of boxes this level= *',ilast-ifirst+1,1)

        if(iffini.eq.1) goto 3000

        ifirst=ilast+1
        ilast=ilastnew

 2000 continue

        ier=1024

 3000 continue

        if(ier.eq.1024) then
        call prinf('adaptive expansion fails with ier=*',ier,1)
        endif

        if(ier.eq.2048) then
        call prinf('number of squares exceed allocated memory *',ier,1)
        endif

        maxlevel=level
        ntot=ilast

        do 3010 i=1,ntot
        c0list(i)=geomlist(1,i)
 3010 continue

        

        call prinf('maxlevel is *',maxlevel,1)
        call prinf('total box number is *',ntot,1)
ccc        call prinf('total root number is *',nrtot,1)
ccc        call prin2_long('all root found *',rootslist,2*nrtot)
ccc        call prin2('root error is *',errlist,nrtot)

c
c remove duplicated roots
c if maxlevel is not 0
c

        if(maxlevel.ne.0) then

         call remove_dup_roots(epsexp,rootslist,nrtot,roots,nr)

        elseif(maxlevel.eq.0) then
         nr=nrtot
         do 3200 i=1,nr
         roots(i)=rootslist(i)
 3200 continue
        endif

        nrtot=0
ccc        call mach_zero(eps)
ccc        eps=sqrt(eps)*10**6
ccc        eps=epsexp*1d10

        do 3210 i=1,nr
        call fun(val,dval,roots(i)) 
ccc        if(abs(val).le.eps) then
ccc        call prin2('abs val is *',abs(val),1)
        rootslist(i)=roots(i)
        nrtot=nrtot+1
        errlist(i)=abs(val)
ccc        endif

 3210 continue


        call prinf('total root number found by the root finder*'
     1          ,nrtot,1)

c
c residue for total root num 
c
        if(ifdval.eq.0) goto 4100

        if(ifres.eq.1) then
         c0=geomlist(1,1)
         sqw=geomlist(2,1)
         call root_number_residue(fun,nres,maxrec,c0,sqw)
         call prinf('total number of roots from residue is *',nres,1)

         ntor=5+log10(sqw)/log10(2*done)

         if(maxrec.gt.ntor) then
           call prinf('Roots  may be near the boundary, 
     1  maxrec is *',maxrec,1)
         endif

        endif


c
c refine with newton
c

        if(ifnewton.eq.1) then
         dzmax=-1
         do 3300 i=1,nrtot
         z0=rootslist(i)

         do 3310 ijk=1,2
         call fun(val,dval,z0)
         dzmag=abs(val/dval)
         if(dzmag.gt.dzmax) dzmax=dzmag
c       call prin2('abs(dz) is *',abs(val/dval),1)
         z0=z0-val/dval
 3310 continue
         rootslist(i)=z0
 3300 continue
         call prin2('max newtown step size is *',dzmax,1)
ccc        call prin2_long('after newton roots*',rootslist,2*nrtot)


c
c update error after newton
c
         do 4000 i=1,nrtot
         call fun(val,dval,rootslist(i)) 
         errlist(i)=abs(val/dval)
 4000 continue
        endif


 4100 continue


        return
        end

        subroutine remove_dup_roots(epsexp,roots1,n1,roots2,n2)
        implicit real *8 (a-h,o-z)
        complex *16 roots1(1),roots2(1)
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c       This subroutine removes duplicated roots between
c       boundaries of neighbouring squares by checking
c       the digits they agree.
c       This approach is relatively robust since
c       multiple roots will be inaccurate.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c       
c
c       input parameters:
c       epsexp:lowerbound of the accuracy of the roots
c       roots1: roots containing redundant ones
c       n1: length of roots1
c
c       output parameters:
c       roots2: roots after redundant ones removed
c       n2: length of roots2
c

        roots2(1)=roots1(1)
        n2=1
        
ccc        call mach_zero(epsexp)
        eps=sqrt(epsexp)*1d-2

        do 2000 i=1,n1
        
        ifput=1
        do 1000 j=1,n2
        d=abs(roots2(j)-roots1(i))
ccc        call prin2('diff is *',d,1)
        if(d.lt.eps) ifput=0
 1000 continue
ccc        call prinf('ifput is *',ifput,1)

        if(ifput.eq.1) then
         n2=n2+1
         roots2(n2)=roots1(i)
        endif

 2000 continue
        
        return
        end



        subroutine divide_square(c0,sqw,c0div,sqwdiv)
        implicit real *8 (a-h,o-z)
        complex *16 c0,c0div(4),ima
        save
        data ima/(0.0d0,1.0d0)/
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c       This subroutine divides a square centered 
c       at c0 with length sqw into four identical  
c       subsquares. Then it returns the centers and 
c       lengths of the four squares.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccc


        
        done=1
        sqwdiv=sqw/(done*2)

        d=sqwdiv/(done*2)
        c0div(1)=c0+(-done+ima)*d
        c0div(2)=c0+(done+ima)*d
        c0div(3)=c0+(-done-ima)*d
        c0div(4)=c0+(done-ima)*d

        return
        end
        

        subroutine complex_root_one_square_adap(epsexp,ifcv,fun,c0,sqw,
     1     vvd,vvoffd,q,z,cw,n,work,croots,norder,roots,nr,ierqr)
        implicit real *8 (a-h,o-z)
        complex *16 vvd(1),vvoffd(1),coefs(10 000),q(1),work(1)
        complex *16 vd(10 000),voffd(10 000)
        complex *16 z(1),cw(1),c0,zz(10 000),cww(10 000)
        complex *16 val,croots(1),roots(1)
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c       This subroutine computes the polynomial expansion coeff.
c       of the input function on a single square, 
c       determining whether the pre-specified accuracy is reached. 
c       If the pre-specified accuracy is not reached, 
c       it returns directly with a flag indicating non-convergence.
c
c       If the accuracy is reached, it proceeds to diagonalize 
c       the corresponding colleague matrix, thus finding all roots 
c       of the function on the square. Then it removes all roots
c       of outside of a slightly extended square of the input one.
c       The remaining roots are returned as the outputs.
cccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       input parameters:
c       epsexp: relative accuracy of the expansion accuracy
c       fun: functions whose roots to be found
c       norder: order of the expansion
c       z: gaussian points on the unit square
c       cw: gaussian weights associated with z
c       n: number of gaussian nodes on each side of the square
c       c0: center of the square
c       sqw: length of the square
c       vvd: diagonal of the colleague matrix
c       vvoffd: off-diagonal of the colleague matrix
c       q: basis matrix (not used)
c       work: stored SVD for least squares
c       No inputs are destroyed by the subroutine
c
c
c       output parameters: 
c       ifcv: if convergent flag
c       croots: all roots of the expansion poly.
c       roots: all roots inside of the extended square
c       nr: length of roots
c
c       ierqr:
c       -1 not convergent, qr not initialted
c       0 qr iteration successful
c       1028 qr fails to converge to acc. set by exiteps
        
c
c
c       work array:
c       vd: copy of vvd
c       voffd: copy of vvoffd
c       zz: scaled/translated gaussian nodes
c       cww: scaled/translated gaussian weights
c
c

cc        ifcv=1 
        done=1
        npt=4*n

        do 100 i=1,norder
        vd(i)=vvd(i)
        voffd(i)=vvoffd(i)
  100 continue


c
c qr not done
c
        ierqr=-1


c
c traslate pts domain into correct domain
c

        do 300 i=1,npt
        zz(i)=z(i)*sqw/(2*done)+c0
        cww(i)=cw(i)*sqw/(2*done)
  300 continue



c
c find expansion coefs
c

        call poly_leastsq_exp(fun,zz,cww,n,q,work,norder,coefs,
     1      sqw)
        
c
c determine the largest |coefs(i)|
c
        coefsmax=-1
        do 450 i=1,norder
        if(coefsmax.lt.abs(coefs(i))) coefsmax=abs(coefs(i))
  450 continue

c
c determine the largest |coefs(i)|
c among the last 11 coefs
c

        coefscut=-1
ccc        istart=norder-10
        istart=norder*4/5
        
ccc        call prinf('start is *',istart,1)
ccc        stop
        do 460 i=istart,norder+1
ccc        coefsave=coefsave+coefs(i)
ccc        call prin2('abs(coef)*',abs(coefs(i)),1)
        if(coefscut.lt.abs(coefs(i))) then 
        coefscut=abs(coefs(i))
        endif

  460 continue


        coefscut=coefscut/coefsmax

c
c   if converge if the relative size
c   of the largest coefs 10 smaller than 1d5*eps
c

        eeps=epsexp**(3*done/4)
ccc        call prin2('eeps is *',eeps,1)
ccc        call prin2('coefscut is *',coefscut,1)

        if(coefscut.gt.eeps) then
ccc          ifcv=3
c         call prinf('norder not enough with norder *',norder,1)
         goto 5000
        endif

        if(ifcv.gt.2) then
         ifcv=ifcv-1
         goto 5000
        endif

         ifcv=ifcv-1

        call sparse_root_struct(croots,vd,voffd,norder,coefs,
     1    ierqr)


        do 500 i=1,norder
        croots(i)=croots(i)*sqw/(2*done)+c0
  500 continue

c        call prinf('tot num of eigenvals found is *',norder,1)

c        call prin2_long('before sort croots are *',croots,norder*2)

c
c estimate max fun value
c

        fmax=0
        do 900 i=1,norder+1
        fmax=fmax+abs(coefs(i))**2
  900 continue
        fmax=sqrt(fmax)
c        call prin2('fmax is *',fmax,1)


        
        nr=0
        do 3000 i=1,norder

ccc        dl=0
ccc        dl=1d-2
ccc        dl=1d-4
        dl=1d-6
ccc        dl=1d-8
        
        ifin=1
        if(abs(real(croots(i)-c0)).ge.(sqw/(2*done)+dl*sqw)) ifin=0
        if(abs(imag(croots(i)-c0)).ge.(sqw/(2*done)+dl*sqw)) ifin=0
        

ccc        if((abs(cy(i)).le.eps).and.(ifin.eq.1)) then
        if(ifin.eq.1) then
cc        if(abs(cy(i)).le.done) then
          nr=nr+1
          roots(nr)=croots(i)
        endif

 3000 continue


 5000 continue

ccc        call prinf('ierqr after one_square is *',ierqr,1)

        return
        end


        subroutine get_pts_and_wts(z,cw,n)
        implicit real *8 (a-h,o-z)
        complex *16 z(1),cw(1)
        dimension ts(10 000),whts(10 000),u(10 000),v(10 000)
        complex *16 s1,s2,s3,s4,ss,cout
        complex *16 c1,c2,c3,c4,c0,c5,c6
        complex *16 ima
        data ima/(0.0d0,1.0d0)/

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      This subroutine constructs Gaussian nodes z and weights cw
c      on the four sides of a square centered at 0 of length 2.
c      Both z and cw are of length 2*4*n, where n is the number of 
c      points on each sides.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        done=1
        

        nview=norder-1
        iflog=0

c cheb nodes and weights on [-1,1]
c
        ifwhts=1
        itype=1
        call legewhts(n,ts,whts,ifwhts)

c
c equally spaced points
c

cccc        h=2/(done*n-1)
cccc
cccc        do 100 i=0,n-1
cccc        ts(i+1)=-1+h*i
ccccccc        whts(i)=h
cccc  100 continue

c
c horizontal length 2*a
c
        a=done
c
c vertical length 2*b
c
        b=done

ccc        call prinf('nc is *',nc,1)
c
c center of the rectangle
c
        c0=0+0*ima
c        c0=1+15*ima

        c1=1*b
        c2=1*a
        c3=1*b
        c4=1*a
ccc        c4=2*a
        c5=2*a
        c6=1*b
        
ccc        call prin2('c5 is *',c5,2)

        call whts_mod(whts,cw,n,c1,c2,c3,c4,c5,c6)
c
        do 1000 i=1,n
c
c C1 contour, 1-i to 1+i
c
        z(i)=b*ts(i)*ima+done*a+c0
c
c C2 contour, 1+i to -1+i
c
        z(n+i)=-a*ts(i)+b*done*ima+c0
c
c C3 contour, -1+i to -1-i
c
        z(2*n+i)=-b*ts(i)*ima-a*done+c0
c
c C4 contour, -1-i to 1-i
c
        z(3*n+i)=a*ts(i)-b*done*ima+c0
c
c C5 contour, -1 to 1
c

ccc        z(4*n+i)=a*ts(i)+c0
c
c C6 contour, -i to i
c

ccc        z(5*n+i)=b*ts(i)*ima+c0

 1000 continue

c
c look at all points
c

        iw=21
        call zquaplot(iw,z,4*n,2,'title*')

        return
        end


        subroutine orth_poly_load(poly,vd,voffd,n)
        implicit real *8 (a-h,o-z)
        real *8 vd(1),voffd(1),poly(1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       This subroutine loads the precomputed polynomial basis
c       with the three-term recurrence coefficients, which
c       form the tridiagonal part of the colleague matrix.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       input parameters:
c       n: order of the polynomail basis
c
c
c       output parameters:
c       vd: diagonal of the tridiagonal matrix, length n
c       voffd: off-diaognal of ......, length n
c       poly: polynomials, length 240*(n+1)

c
c n is the order of expansion
c
        call diag_elem(vd,n) 
        call subdiag_elem(voffd,n) 
c
c n order polynomials have n+1 terms
c
        call orth_vecs(poly,n+1) 

        return
        end
        
        
        
        subroutine whts_mod(whts,w,n,c1,c2,c3,c4,c5,c6)
        implicit real *8 (a-h,o-z)
        dimension whts(1)
        complex *16 w(1)
        complex *16 ima,c2,c3,c4,c5,c6
        save
        data ima/(0.0d0,1.0d0)/

        done=1
        j=0
ccc        c1=done
        do 900 i=1,n
        w(i+j*n)=whts(i)*c1
  900 continue

        j=1
ccc        c2=done
        do 1000 i=1,n
        w(i+j*n)=whts(i)*c2
 1000 continue

        j=2
ccc        c3=done
        do 1100 i=1,n
        w(i+j*n)=whts(i)*c3
 1100 continue

        j=3
ccc        c4=done*3
        do 1200 i=1,n
        w(i+j*n)=whts(i)*c4
 1200 continue
        
cccc        j=4
cccc        do 1300 i=1,n
cccc        w(i+j*n)=whts(i)*c5
cccc 1300 continue
cccc
cccc        j=5
cccc        do 1400 i=1,n
cccc        w(i+j*n)=whts(i)*c6
cccc 1400 continue

        return
        end




        subroutine poly_leastsq_exp(fun,z,cw,n,q,w,norder,coeff,
     1     sqw)
        implicit real *8 (a-h,o-z)
        complex *16 q(4*n,norder+1),coeff(1),a(4*n,norder+1)
        complex *16 z(1),cw(1),zfg(1)
        complex *16 cval(1000),val,fval(1000)
        complex *16 w(1)
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/

        m=norder+1
        npt=4*n
        done=1
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c This subroutine computes expansion coefficients of the input fun
c in the precomputed basis q. It does so by firstly adding weights
c sqrt(cw_i) to values of fun at points z_i, forming the cval vector,
c and adding sqrt(cw_i) to each basis at z_i, forming the a matrix.
c The coefficients are computed by solving the least squares problem
c               a*coeff  =  cval
c via the SVD approach.
c
c Array w contains the precomputed SVD of a, so ncleast2 can
c be called directly for computing coeff.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c
c solve coeff for min ||cval - a*coeff||
c

ccc        fmax=0

        do 1000 i=1,npt
ccc        call test_cfun1(fval(i),z(i))
        call fun(fval(i),val,z(i))
ccc        if(fmax.lt.abs(fval(i))) fmax=abs(fval(i))
        cval(i)=sqrt(cw(i))*fval(i)
ccc        cval(i)=fval(i)
 1000 continue

ccc        do 1200 j=1,m
cccc        do 1100 i=1,npt
ccc        do 1100 i=1,4*n
ccc        a(i,j)=q(i,j)*sqrt(cw(i))
cccccc        a(i,j)=q(i,j)
ccc 1100 continue
ccc 1200 continue
ccc        call mach_zero(eps)
ccc        eps=eps*5
ccc        call ncleastsq(a,npt,m,eps,ncols,rnorms,w) 
cccccc        call prinf('ncol is *',ncols,1)
ccc
        call ncleasts2(w,cval,coeff)



        return
        end




cccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccc
c
c eigenvale solver start here
c
cccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccc

        subroutine sparse_root_struct(croot,vd,voffd,norder,coefs,
     1   ierqr)
        implicit real *8 (a-h,o-z)
        complex *16 croot(norder)
        complex *16 vd(1),voffd(1),coefs(1),zfg(1)
        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   This subroutine computes the roots of a polynomial 
c                     norder
c           p(z)= \Sum       coefs_(i+1) P_i            (1)
c                     i=0
c   with a user specified coefs.
c
c   The vectors vd, voffd specify the three-term recurrence  that define 
c   polynomials {P_i}_{i=0}^norder. 
c
c   It does so by finding the eigenvalues of the so-called 
c   colleague matrix, constructed from the expansion coefficients. 
c   This subroutine does not assume {P_i} to classical polynomials;
c   it only assume they satisfy a three-term recurrence, so
c   the colleague matrix is in general complex symmetric, as opposed 
c   to being real symmetric. As a result, complex orthogonal transforms
c   are used to diagonalize the colleague matrix.   
c
c   The claim to fame of this subroutine is that it has cost O(n^2).
c   Besides, numerical experiments strongly suggest that it is 
c   also componenetwise backward stable in the sense that the computed
c   eigenvalues are the exact eigenvalues of a perturbed matrix, where    
c   the perturbation is proportional to the size of each entry * eps.
c
c
c                          WARNINGS:
c   1. Detailed perturbation analysis has not been done.
c   2. This subroutine uses complex orthogonal transforms for 
c      diagonalizing the colleague matrix, so the numerical stability
c      is not guaranteed!
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                               inputs:
c   norder: order of the polynomial expansion
c   vd: diagonal of the colleague matrix, contains norder elements,
c       not touched by the subroutine
c   voffd: superdiagonal of the colleague matrix, contains norder 
c          elements not touched by the subroutine
c   coefs: polynomial expansion coefficients, contains norder elements,
c          not touched by the subroutine
c
c                               outputs:
c   croot: computed eigenvalues of the colleague matrix, thus
c          the roots of the polynomial expansion (1), contains
c          norder elements
c   ierqr:
c        0 qr iteration successful
c        1028 qr fails to converge to acc. set by exiteps
        
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        
        done=1
        
ccc        call prin2('in wilk vd is *',vd,2*norder)
ccc        call prin2('in wilk voffd is *',voffd,2*norder)
ccc        call prin2('in wilk coefs is *',coefs,2*norder)

        call mach_zero(exiteps)
cccc        exiteps=exiteps*10

        csize=0
        do 1000,i=1,norder-1
        csize=abs(vd(i))**2 + 2*abs(voffd(i))**2 + csize
 1000 continue

        csize=csize+abs(vd(norder))**2

        csize=sqrt(csize)

cc        call prin2('csize is *',csize,1)

        exiteps=exiteps*csize

ccc        call prin2('exit eps in wilk is *',exiteps,1)

        call cqr_eigenvals_wilk(croot,norder,exiteps,ierqr,
     1              vd,voffd,coefs)


        return
        end



        subroutine cqr_one_iter_dumb(n,vd,voffd,pv,qv)
        implicit real *8 (a-h,o-z)
        complex *16 vd(1),voffd(1),qv(1),coefs(1),pv(1)
        complex *16 cpqv(10 00),gv(10 00)
        complex *16 sprot(2,n),c,s,cc,cs
        complex *16 d0,d1,d2,c0,s0
ccc        dimension rotnorm(10 000),aaa(10 000)
c
c coefs is not used in the old code
c


c
c sparse qr starts here
c
        do 10 i=1,n
        cpqv(i)=qv(i)
        gv(i)=voffd(i)
   10 continue

        do 1300 k=n,2,-1

        d1=voffd(k-1)+pv(k-1)*qv(k)
        d2=vd(k)+pv(k)*qv(k)

ccc        call prin2('rot norm is *',rotnorm(k-1),1)
        call find_rot_dumb(d1,d2,c,s)
        sprot(1,k-1)=c
        sprot(2,k-1)=s



c
c apply sparse row rotation
c 

        c=sprot(1,k-1)
        s=sprot(2,k-1)

        

        cc=c
        cs=s


c
c sparse rotate sub and subsub
c
        if(k.gt.2) then
          d1=c*gv(k-2)-s*(-cpqv(k)*pv(k-2))

          gv(k-2)=d1
        endif
c
c sparse rotate sub and diag
c
        d1=c*vd(k-1)-s*gv(k-1)
        d2=cs*vd(k-1)+cc*gv(k-1)

        vd(k-1)=d1
        gv(k-1)=d2

c
c sparse rotate diag and sup
c
        d1=c*voffd(k-1)-s*vd(k)
        d2=cs*voffd(k-1)+cc*vd(k)


        voffd(k-1)=d1
        vd(k)=d2
c
c rotate pv
c
        d1=c*pv(k-1)-s*pv(k)
        d2=cs*pv(k-1)+cc*pv(k)
        
        pv(k-1)=d1
        pv(k)=d2
        
c
c apply correction
c
        
        ddd1=abs(pv(k-1)*qv(k))**2+abs(pv(k)*qv(k))**2
        ddd2=abs(voffd(k-1))**2+abs(vd(k))**2
        if(ddd1.ge.ddd2) then
            pv(k-1)=-voffd(k-1)/qv(k)
        endif

        
c
c rotate cpqv
c
        d1=c*cpqv(k-1)-s*cpqv(k)
        d2=cs*cpqv(k-1)+cc*cpqv(k)
        cpqv(k-1)=d1
        cpqv(k)=d2


 1300 continue

        
c
c sparse rotate column
c
        do 2200 k=n,2,-1



        c=sprot(1,k-1)
        s=sprot(2,k-1)
        cc=c
        cs=s

c
c rotate diag and sup column
c
        d1=cc*vd(k-1)-cs*(-pv(k-1)*qv(k))
        d2=s*vd(k-1)+c*(-pv(k-1)*qv(k))

        vd(k-1)=d1
        voffd(k-1)=d2


c
c rotate sub and diag column
c
        d1=cc*gv(k-1)-cs*vd(k)
        d2=s*gv(k-1)+c*vd(k)
        gv(k-1)=d1
        vd(k)=d2

c
c rotate qv
c

        d1=c*qv(k-1)-s*qv(k)
        d2=cs*qv(k-1)+cc*qv(k)

        qv(k-1)=d1
        qv(k)=d2


 2200 continue

        return
        end


        subroutine find_rot_dumb(x1,x2,c,s)
        implicit real *8 (a-h,o-z)
        complex *16 x1,x2,c,s,v,vv,ss,cc,d
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c   This subroutine finds the 2*2 complex orthogonal transform Q,
c   so that
c           [c -s] [x1]  =  [      0          ]
c           [s  c] [x2]     [sqrt(x1**2+x2**2)]
c
c                       WARNING:
c   Size of c and s are not bounded! They can be large.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        done=1
        pi=atan(done)*4
        v=sqrt((x1)**2+(x2)**2)
        vv=sqrt(abs(x1)**2+abs(x2)**2)
        if(abs(vv).eq.0) then
         c=1
         s=0
         goto 1000 
        endif

        c=x2/v 

        s=x1/v

 1000 continue

c
c if |x1|**2+|x2|**2 is not small
c

c        call prin2('c is *',c,2)
c        call prin2('s is *',s,2)

        return
        end
        

        subroutine eigenvals_quadr_sol2_elem(a11,a21,a12,a22,
     1    clam1,clam2)
        implicit real *8 (a-h,o-z)
        save
        complex *16 aa(2,2),clam1,clam2,discr,a11,a21,a12,a22
        data four/4.0d0/,half/0.5d0/

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       This subroutine computes the two eigenvalues of a complex
c       2 by 2 matrix whose elements are a11,a21,a12,a22. 
c
c       clam1 and clam2 are the two eigenvalues.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
        aa(1,1)=a11
        aa(1,2)=a12
        aa(2,1)=a21
        aa(2,2)=a22
        
        
        discr=sqrt( (aa(1,1)+aa(2,2))**2+
     1      four*( aa(1,2)*aa(2,1) -aa(1,1)*aa(2,2)) )
c
        clam1=aa(1,1)+aa(2,2)+discr
     1      
        clam1=clam1*half
c
        clam2=aa(1,1)+aa(2,2)-discr
        clam2=clam2*half
c
        return
        end


        subroutine cqr_eigenvals_wilk(clam,n,exiteps,ier,
     1              vd,voffd,coefs)
        implicit real *8 (a-h,o-z)
        dimension nstack(10000,3)
        complex *16 clam(1),shift,coefs(1)
        complex *16 vd(1),voffd(1),pv(10 00),qv(10 00),clam1,clam2
        complex *16 a11,a12,a21,a22,clam3,offelem,elemmin
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       This subroutine computes the n roots of a polynomial of order n 
c       whose coefficients are contained in the coefs vector. It is done
c       by finds the spectrum of the generalized colleague matrix of the 
c       form 
c
c                    C = A + pv \circ qv^T,                            (1)
c
c       where A is a complex symmetric matrix, and pv, qv are a pair of vectors.
c       It is done by a special complex orthogonal QR with wilkinson
c       shifts.
c
c       The claim to fame of this subroutine is that it costs O(n^2)
c       operations.
c       Numerical experiments suggest that it has a special form 
c       of backward stability, where
c       the computed eigenvalues are the exact eigenvalues of the matrix 
c
c                   (A+dA) + (pv+dpv) \circ (qv+dqv)^T                  (2)
c
c       with dA \approx eps||A||, dpv \approx eps||pv||, and 
c       dqv \approx eps||qv||.
c
c       PLEASE NOTE that this is a very specialized subroutine, in 
c       in that the matrix C is BOTH of the form (1) and lower 
c       Hessenberg, so that 
c
c                    C(i,j)=0                                           (3)
c
c       for all 
c
c                    j > i+1.                                           (4)
c
c       Note also that, for matrices of this form, A is determined
c       entirely by its diagonal and superdiagonal together with the
c       vectors pv and qv.
c       
c                               WARNING:
c       Unlike usual QR algorithms, complex orthogonal rotations
c       are used, instead of unitary ones, in order to maintain
c       the complex symmetry of A. As a result, numerical stability
c       is not guaranteed and breakdowns can occur. However, breakdowns
c       have not been encountered once shifts are introduced; the QR
c       converges fast enough before breakdowns happen.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c                      input parameters:
c
c   vd: diagonal of A, contains n elements, destroyed by the
c       subroutine.
c   voffd: contains n elements, voffd(1),...,voffd(n-1) are the  
c          superdiaognal of A, and voffd(n) is used to form qv with coefs,
c          destroyed by the subroutine.
c   coefs: (n+1) coefficients of a polynomial of order n, used to 
c          form the q vector, q= coefs(1:n)/(coefs(n+1)*voffd(n))
c   exiteps: the ABSOLUTE precision to which the eigenvalues will be
c       determined. 
c
c
c
c                      output parameters:
c   clam:  eigenvalues of the matrix B of n complex numbers
c   ier:
c        0 successful
c        1028 qr fails to converge to acc. set by exiteps
c        
c
c                      working arrays:
c   pv: the vector pv, must be at least n complex numbers
c   qv: the vector qv, must be at least n complex numbers
c   nstack: stack variable for splitting in QR, depth at most 10000
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


ccc        call prinf('n is cqr *',n,1)
        do 100 i=1,n
        qv(i)=voffd(n)*coefs(i)/coefs(n+1)
        pv(i)=0
  100 continue

        pv(n)=-1

ccc        call prin2('in wilkqr pv is *',pv,2*n)
ccc        call prin2('in wilkqr qv is *',qv,2*n)

        iternum=0
        is=1 
        maxstack=is
        nstack(is,1)=1
        nstack(is,2)=n
        nstack(is,3)=nstack(is,2)-nstack(is,1)+1
cccc        call prinf('current size is *',nstack(1,3),1)

        
        ier=0
c        call prini(6,14)
        do 5000 ijk=1,10000
        if(is.eq.0) goto 6000

        nn=nstack(is,3)
        if((nn.eq.1).or.(nn.eq.2)) goto 2000

        do 1000 i=nstack(is,1),(nstack(is,2)-1)
        offelem=voffd(i)+pv(i)*qv(i+1)

        if(abs(offelem).le.exiteps) then
c        call prin2('supdiag elem is *',offelem,2)
c        call prinf('at position *',i,1)

c
c go down the stack
c 
        is=is+1
        nstack(is,1)=nstack(is-1,1)
        nstack(is,2)=i
        nstack(is,3)=nstack(is,2)-nstack(is,1)+1

        nstack(is-1,1)=i+1
        nstack(is-1,3)=nstack(is-1,2)-nstack(is-1,1)+1

        if(maxstack.lt.is) maxstack=is

c        call prinf('current stack from *',nstack(is,1),1)
c        call prinf('to *',nstack(is,2),1)
c        call prinf('of size *',nstack(is,3),1)

c        call prinf('previous stack from *',nstack(is-1,1),1)
c        call prinf('to *',nstack(is-1,2),1)
c        call prinf('of size *',nstack(is-1,3),1)
        goto 2000   

        endif
        
 1000 continue
ccc        call prin2('min supdiag elem is *',elemmin,2)
ccc        call prinf('at position *',ip,1)
 2000 continue
c        call prinf('stack number is *',is,1)

        ii=nstack(is,1)
        nn=nstack(is,3)
        
        if(nn.eq.1) then
        a11=vd(ii)+pv(ii)*qv(ii)
        clam(ii)=a11
c        call prin2('size 1, clam is *',a11,2)

c
c go up the stack
c 
        is=is-1

        goto 3000
        endif

        if(nn.eq.2) then

        a11=vd(ii)+pv(ii)*qv(ii)
        a12=voffd(ii)+pv(ii)*qv(ii+1)
        a21=voffd(ii)+pv(ii+1)*qv(ii)
        a22=vd(ii+1)+pv(ii+1)*qv(ii+1)

        call eigenvals_quadr_sol2_elem(a11,a21,a12,a22,
     1    clam1,clam2)
        clam(ii)=clam1
        clam(ii+1)=clam2

c        call prin2('size 2, clam1 is *',clam1,2)
c        call prin2('clam2 is *',clam2,2)

c
c go up the stack
c 
        is=is-1

        goto 3000
        endif

        call qr_one_step(nn,vd(ii),voffd(ii),pv(ii),qv(ii))
        iternum=iternum+1

 3000 continue

 5000 continue
        ier=1028
 6000 continue

c        call prini(6,14)
ccc        call prin2('clam is *',clam,2*n)
c        call prinf('cqr wilk ier is *',ier,1)
c        call prinf('max stack level is *',maxstack,1)
c        call prinf('cqr wilk iteration num is *',iternum,1)
c        call prini(6,13)

cc        if(ier.ne.0) then 
cc        call prinf('cqr wilk fails with ier=*',ier,1)
cc        endif


        if(ier.eq.1028) then 
        call prinf('cqr wilk fails to converge 
     1            after 1000 iterations, ier=* ',ier,1)
        endif




        return
        end

        subroutine qr_one_step(n,vd,voffd,pv,qv)
        implicit real *8 (a-h,o-z)
        complex *16 vd(1),voffd(1),pv(1),qv(1)
        complex *16 clam1,clam2,shift,a11,a12,a21,a22

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c       This subroutine performs one (complex orthogonal) QR iteration
c       to the matrix represented by generators vd, voffd, pv, qv, 
c       all containing n complex numbers. The generators are updated
c       by the subroutine.
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


        a11=vd(1)+pv(1)*qv(1)
        a12=voffd(1)+pv(1)*qv(2)
        a21=voffd(1)+pv(2)*qv(1)
        a22=vd(2)+pv(2)*qv(2)
ccc eigenvals_quadr_sol2_elem(a11,a21,a12,a22,
ccc     1    clam1,clam2)
        call eigenvals_quadr_sol2_elem(a11,a21,a12,a22,
     1    clam1,clam2)

        shift=clam1
        if(abs(clam1-a11) .gt. abs(clam2-a11) ) shift=clam2


ccc        call prin2('clam1 is *',clam1,2)
ccc        call prin2('clam2 is *',clam2,2)
ccc        call prin2('a11 is *',a11,2)
ccc        call prin2('shift is *',shift,2)

c
c apply shift
c
        do 1000 i=1,n
        vd(i)=vd(i)-shift
 1000 continue

        call cqr_one_iter_dumb(n,vd,voffd,pv,qv)

        do 1100 i=1,n
        vd(i)=vd(i)+shift
 1100 continue


        return
        end



        subroutine root_number_residue(fun,nres,mrec,c0,sqw)
        implicit real *8 (a-h,o-z)
        complex *16 rint1,rint2,rint3,rint4,rint,val
        complex *16 trap1,trap2,trap3,trap4
        complex *16 c0,par(2)
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/
ccc        external fun
        external test_cfun_int1
        external test_cfun_int2
        external test_cfun_int3
        external test_cfun_int4
        
c
c residue integration
c to find num of roots
c
        done=1
ccc        call prini(6,19)
        b=sqw/(2*done)
        a=-b
cc        call prin2('b is *',b,1)
cc        call prin2('a is *',a,1)
cc        call prin2('center c0 is *',c0,2)
        
        par(1)=c0
        par(2)=sqw
ccc        a=-1
ccc        b=1

        call mach_zero(eps)
        eps=sqrt(eps)
ccc        eps=1.0d-15
        m=32
        np=1

        mrec=0
        call cadapgau_new(ier,a,b,test_cfun_int1,fun,par,m,eps,
     1      rint1,maxrec,numint)
ccc        call prinf('ier is *',ier,1)
ccc        call prin2('rint1 is *',rint1,2)
        if(ier.ne.0) then
        call prinf('C1 int failure, ier is *',ier,1)
        endif
ccc        call prinf('C1 maxrec is *',maxrec,1)
        if(mrec.lt.maxrec) mrec=maxrec

        np=2
        call cadapgau_new(ier,a,b,test_cfun_int2,fun,par,m,eps,
     1      rint2,maxrec,numint)
ccc        call prinf('ier is *',ier,1)
ccc        call prin2('rint2 is *',rint2,2)

        if(ier.ne.0) then
        call prinf('C2 int failure, ier is *',ier,1)
        endif
ccc        call prinf('C2 maxrec is *',maxrec,1)
        if(mrec.lt.maxrec) mrec=maxrec

        np=3
        call cadapgau_new(ier,a,b,test_cfun_int3,fun,par,m,eps,
     1      rint3,maxrec,numint)
ccc        call prinf('ier is *',ier,1)
ccc        call prin2('rint3 is *',rint3,2)

        if(ier.ne.0) then
        call prinf('C3 int failure, ier is *',ier,1)
        endif
ccc        call prinf('C3 maxrec is *',maxrec,1)
        if(mrec.lt.maxrec) mrec=maxrec

        np=4
        call cadapgau_new(ier,a,b,test_cfun_int4,fun,par,m,eps,
     1      rint4,maxrec,numint)
ccc        call prinf('ier is *',ier,1)
ccc        call prin2('rint4 is *',rint4,2)

        if(ier.ne.0) then
        call prinf('C4 int failure, ier is *',ier,1)
        endif
ccc        call prinf('C4 maxrec is *',maxrec,1)
        if(mrec.lt.maxrec) mrec=maxrec

        

        rint=rint1*ima-rint2-rint3*ima+rint4
ccc        call prin2('contour int val is *',rint,2)
        pi=4*atan(done)
ccc        call prin2_long('num of roots from residue is *',rint/(2*pi),2)

        nres=imag(rint/(2*pi))+0.1
cccc        if(abs(nres*done-imag(rint/(2*pi))).ge.0.5) nres=nres+1

ccc        call prinf('num of roots integer is *',nres,1)

ccc        call prini(6,13)

ccc        call prinf('maxrec is *',mrec,1)

        return
        end



        subroutine copy_cmat(a,b,n)
        implicit real *8 (a-h,o-z)
        complex *16 a(n,n),b(n,n)
        
        do 2000 i=1,n
        do 1000 j=1,n
        b(i,j)=a(i,j)
 1000 continue
 2000 continue
        
        return
        end


        subroutine test_cfun_int1(x,fun,par,val)
        implicit real *8 (a-h,o-z)
        complex *16 val,val1,val2,z
        complex *16 par(1),c0
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/


        done=1

        sqw=real(par(2))/(2*done) 
        c0=par(1)
        z=(done*sqw+ima*x)+c0

        call fun(val1,val2,z)
ccc        call dfun(val2,z)

        val=val2/val1

        return
        end

        subroutine test_cfun_int2(x,fun,par,val)
        implicit real *8 (a-h,o-z)
        complex *16 val,val1,val2,z
        complex *16 par(1),c0
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/


        done=1

        sqw=real(par(2))/(2*done) 
        c0=par(1)
        z=(done*x+ima*sqw)+c0

        call fun(val1,val2,z)

ccc        call fun(val1,z)
ccc        call dfun(val2,z)

        val=val2/val1

        return
        end
        
        subroutine test_cfun_int3(x,fun,par,val)

        implicit real *8 (a-h,o-z)
        complex *16 val,val1,val2,z
        complex *16 par(1),c0
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/


        done=1

        sqw=real(par(2))/(2*done) 
        c0=par(1)
        z=(-done*sqw+ima*x)+c0

        call fun(val1,val2,z)

ccc        call fun(val1,z)
ccc        call dfun(val2,z)

        val=val2/val1

        return
        end

        subroutine test_cfun_int4(x,fun,par,val)
        implicit real *8 (a-h,o-z)
        complex *16 val,val1,val2,z
        complex *16 par(1),c0
        complex *16 ima
        save
        data ima/(0.0d0,1.0d0)/


        done=1

        sqw=real(par(2))/(2*done) 
        c0=par(1)
        z=(done*x-ima*sqw)+c0

        call fun(val1,val2,z)
ccc        call fun(val1,z)
ccc        call dfun(val2,z)

        val=val2/val1

        return
        end



