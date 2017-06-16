! ==============================================================================
! storm_module
!
! Module contains routines for constructing a wind and pressure field based on
! the Holland-Silva hurricane module.
!
! Many of these routines are based loosely on PADCIRC version 45.12 03/17/2006
! ==============================================================================
module silva_storm_module

    implicit none
    save

    ! silva storm type definition
    type silva_storm_type
        ! Fore/hindcast size and current position
        integer :: num_casts

        ! These parameters are located at time points but are interpolated in
        ! time and space when the relevant fields are requested.

        ! Location of storm
        ! Track is a triplet with (time,longitude,latitude)
        real(kind=8), allocatable :: track(:,:)

        ! Storm parameterization
        real(kind=8), allocatable :: max_wind_radius(:)
        real(kind=8), allocatable :: max_wind_speed(:)
        real(kind=8), allocatable :: central_pressure(:)
        real(kind=8), allocatable :: rrp(:)

        ! Approximate velocity of storm, approximated via the track points
        ! using a first order difference on the sphere
        real(kind=8), allocatable :: velocity(:,:)

    end type silva_storm_type

    logical, private :: module_setup = .false.

    ! Interal tracking variables for storm
    real(kind=8), private :: A,B
    integer, private :: last_storm_index

    logical, private :: DEBUG = .false.

    ! Atmospheric boundary layer, input variable in ADCIRC but always is
    ! set to the following value
    real(kind=8), parameter :: atmos_boundary_layer = 0.9d0

    ! Sampling adjustment from 1 min to 10 min winds
    real(kind=8), parameter :: sampling_time = 0.88d0

    ! Storm field ramping width - Represents crudely the ramping radial area
    real(kind=8), parameter :: RAMP_WIDTH = 100.0d3

contains

    ! Setup routine for the silva model
    subroutine set_silva_storm(storm_data_path, storm, log_unit)

        use geoclaw_module, only: deg2rad, spherical_distance, coordinate_system
        use amr_module, only: t0, rinfinity

        implicit none

        ! Subroutine I/O
        character(len=*), optional :: storm_data_path
        type(silva_storm_type), intent(in out) :: storm
        integer, intent(in) :: log_unit

        ! Local storage
        integer, parameter :: data_file = 701
        integer :: i, k, io_status, num_casts
        real(kind=8) :: forecast_time,last_time,x(2),y(2),ds,dt,dx,dy,theta

        ! Reading buffer variables
        integer :: year,month,day,hour,forecast,lat,lon,max_wind_speed
        integer :: central_pressure,RRP,max_wind_radius
        character(len=4) :: cast_type
        character(len=1) :: direction(2)

        ! Note that the JAM file format has not been tested yet and will later
        ! be added as an option.
        character(len=4), parameter :: file_format = "NOAA"

        ! File format string
        character(len=*), parameter :: JMA_FORMAT = "(i2,i2,i2,i2,8x,i3,1x,"//&
                            "i4,1x,i4,5x,i3)"
        character(len=*), parameter :: NOAA_FORMAT = "(8x,i4,i2,i2,i2,6x,a4,"//&
                            "2x,i3,1x,i4,a1,2x,i4,a1,2x,i3,2x,i4,47x,i3,2x,i3)"

        if (.not. module_setup) then

            ! Storm type only works on lat-long coordinate systems
            if (coordinate_system /= 2) then
                stop "silva storm type does only works on lat-long coordinates."
            endif

            ! Open data file
            if (present(storm_data_path)) then
                print *,'Reading storm date file ',storm_data_path
                open(unit=data_file,file=storm_data_path,status='old', &
                     action='read',iostat=io_status)
            else
                print *,'Reading storm date file ./storm.data'
                open(unit=data_file,file="./storm.data",status='old', &
                     action='read',iostat=io_status)
            endif
            if (io_status /= 0) then
                print "(a,i2)", "Error opening storm data file. status = ", io_status
                stop
            endif

            ! Count number of data lines
            num_casts = 0
            last_time = -rinfinity
            do
                if (file_format == "NOAA") then
                    read(data_file,fmt=NOAA_FORMAT,iostat=io_status) year,month,day, &
                        hour,cast_type,forecast,lat,direction(2),lon,direction(1), &
                        max_wind_speed,central_pressure,RRP,max_wind_radius
                else if (file_format == "JAM") then
                    ! JAM may be missing RRP parameter, may need to set this based
                    ! on other data in the file.  It is only used in the field
                    ! ramping function so it might not be an issue
                    read(data_file,fmt=JMA_FORMAT,iostat=io_status) year, month, day, &
                            hour, lat, lon, central_pressure, max_wind_speed, max_wind_radius
                else
                    print *,"ERROR - Unrecognized storm data file format."
                    stop
                endif

                ! Exit loop if we ran into an error or we reached the end of the file
                if (io_status /= 0) exit

                ! Skip counting this line if time is repeated
                forecast_time = date_to_seconds(year,month,day,hour,0,0.d0)
                if (abs(forecast_time - last_time) >= 1.8d3) then
                    num_casts = num_casts + 1
                endif
                last_time = forecast_time
            end do
            rewind(data_file)

            write(log_unit,"('Forecasts = ',i3)") num_casts

            ! Allocate storm parameter file variables
            allocate(storm%track(3,num_casts))
            allocate(storm%max_wind_speed(num_casts))
            allocate(storm%max_wind_radius(num_casts))
            allocate(storm%central_pressure(num_casts))
            allocate(storm%rrp(num_casts))

            ! Now re-read the file's contents
            i = 0
            do while (i < num_casts)
                if (file_format == "NOAA") then
                    read(data_file,fmt=NOAA_FORMAT) year,month,day,hour,cast_type, &
                        forecast,lat,direction(2),lon,direction(1),max_wind_speed, &
                        central_pressure,RRP,max_wind_radius
                else if (file_format == "JAM") then
                    read(data_file,fmt=JMA_FORMAT,iostat=io_status) year, month, day, &
                            hour, lat, lon, central_pressure, max_wind_speed, max_wind_radius
                else
                    print *,"ERROR - Unrecognized storm data file format."
                    stop
                endif


                ! Skip counting this line if time is repeated
                forecast_time = date_to_seconds(year,month,day,hour,0,0.d0)
                if (abs(forecast_time - last_time) < 1.8d3) then
                    cycle
                endif
                i = i + 1
                last_time = forecast_time

                ! Storm position
                ! Conversions:
                !  lon - Convert 10ths of degs to degs, depends on E,W
                !  lat - Convert 10ths of degs to degs
                storm%track(1,i) = date_to_seconds(year,month,day,hour,0,0.d0)
                if (direction(1) == "E") then
                    storm%track(2,i) = real(lon,kind=8) / 10.d0
                else
                    storm%track(2,i) = -real(lon,kind=8) / 10.d0
                endif
                if (direction(2) == "N") then
                    storm%track(3,i) = real(lat,kind=8) / 10.d0
                else
                    storm%track(3,i) = -real(lat,kind=8) / 10.d0
                endif

                ! Storm intensity
                ! Conversions:
                !  max_wind_speed - Convert knots to m/s
                !  NO: max_wind_radius  - convert from nm to m
                !  max_wind_radius  - cyclostrophic radius (Silva) [km] to m
                !  NO: central_pressure - convert from mbar to Pa
                !  Radius of last isobar contour - convert from nm to m
                storm%max_wind_speed(i) = real(max_wind_speed,kind=8) * 0.51444444d0
                storm%max_wind_radius(i) = (0.4785* real(central_pressure,kind=8) - 413.01)*1000.0
                storm%central_pressure(i) = real(central_pressure,kind=8) * 100.d0
                storm%rrp(i) = real(RRP,kind=8) * 1.852000003180799d0

            enddo

            ! Calculate storm speed
            allocate(storm%velocity(2,num_casts))
            do i=1,num_casts - 1
                ! Calculate velocity based on great circle distance between

                ! locations of storm
                x = storm%track(2:3,i)
                y = storm%track(2:3,i+1)

                dt = storm%track(1,i + 1) - storm%track(1,i)

                ds = spherical_distance(x(1), 0.5d0 * (x(2) + y(2)), &
                                        y(1), 0.5d0 * (x(2) + y(2)))
                storm%velocity(1,i) = sign(ds / dt,y(1) - x(1))


                ds = spherical_distance(0.5d0 * (x(1) + y(1)), x(2), &
                                        0.5d0 * (x(1) + y(1)), y(2))
                storm%velocity(2,i) = sign(ds / dt,y(2) - x(2))
            end do

            ! Use last approximation for velocity point going forward
            storm%velocity(:,num_casts) = storm%velocity(:,num_casts - 1)

            ! Record number of casts
            storm%num_casts = num_casts

            if (t0 < storm%track(1,1)) then
                print *,t0,storm%track(1,1)
                stop "Start time is before first forecast time."
            endif

            ! This is used to speed up searching for correct storm data
            last_storm_index = 2
            last_storm_index = storm_index(t0,storm)
            if (last_storm_index == -1) then
                print *,"Forecast not found for time ",t0,'.'
                stop
            endif

            ! Log everything to the surge log file
            write(log_unit,*) ""
            write(log_unit,*) "Storm Track and Strength"
            write(log_unit,*) ""
            do i=1,storm%num_casts
                write(log_unit,"(8e26.16)") (storm%track(k,i),k=1,3),  &
                                            (storm%velocity(k,i),k=1,2), &
                                             storm%max_wind_radius(i), &
                                             storm%max_wind_speed(i),  &
                                             storm%central_pressure(i)
            enddo

            module_setup = .true.
        end if

    end subroutine set_silva_storm


    ! ==========================================================================
    !  real(kind=8) pure date_to_seconds(year,months,days,hours,minutes,seconds)
    !    Convert time from year, month, day, hour, min, sec to seconds since the
    !    beginning of the year.
    ! ==========================================================================
    pure real(kind=8) function date_to_seconds(year,months,days,hours,minutes, &
                                               seconds) result(time)

        implicit none

        ! Input
        integer, intent(in) :: year, months, days, hours, minutes
        real(kind=8), intent(in) :: seconds

        ! Local storage
        integer :: total_days

        ! Count number of days
        total_days = days

        ! Add days for months that have already passed
        if (months > 1) total_days = total_days + 31
        if (months > 2) then
            if (int(year / 4) * 4 == year) then
                total_days = total_days + 29
            else
                total_days = total_days + 28
            endif
        endif
        if (months > 3)  total_days = total_days + 31
        if (months > 4)  total_days = total_days + 30
        if (months > 5)  total_days = total_days + 31
        if (months > 6)  total_days = total_days + 30
        if (months > 7)  total_days = total_days + 31
        if (months > 8)  total_days = total_days + 31
        if (months > 9)  total_days = total_days + 30
        if (months > 10) total_days = total_days + 31
        if (months > 11) total_days = total_days + 30

        ! Convert everything to seconds since the beginning of the year
        time = real((total_days - 1) * 86400 + hours * 3600 + minutes * 60,kind=8)
        time = time + seconds

    end function date_to_seconds


    ! ==========================================================================
    !  silva_storm_location(t,storm)
    !    Interpolate location of hurricane in the current time interval
    ! ==========================================================================
    function silva_storm_location(t,storm) result(location)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(silva_storm_type), intent(in out) :: storm

        ! Output
        real(kind=8) :: location(2)

        ! Junk storage
        real(kind=8) :: junk(2)

        call get_silva_storm_data(t,storm,location, &
                                        junk,junk(1),junk(1),junk(1),junk(1))

    end function silva_storm_location

    ! ==========================================================================
    !  silva_storm_direction
    !   Angle off of due north that the storm is traveling
    ! ==========================================================================
    real(kind=8) function silva_storm_direction(t, storm) result(theta)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(silva_storm_type), intent(in) :: storm

        ! Locals
        real(kind=8) :: junk(2), velocity(2)

        ! Fetch velocity of storm which has direction encoded in it
        call get_silva_storm_data(t, storm, junk, velocity, junk(1),  &
                                                    junk(1), junk(1), junk(1))

        ! Unit directional vector
        theta = atan2(velocity(2),velocity(1))

    end function silva_storm_direction

    ! ==========================================================================
    !  storm_index(t,storm)
    !    Finds the index of the next storm data point
    ! ==========================================================================
    integer pure function storm_index(t,storm) result(index)

        implicit none

        ! Input
        real(kind=8), intent(in) :: t
        type(silva_storm_type), intent(in) :: storm

        ! Locals
        real(kind=8) :: t0,t1
        logical :: found

        ! Figure out where we are relative to the last time we checked for the
        ! index (stored in last_storm_index)

        ! Check if we are already beyond the end of the last forecast time
        if (last_storm_index == storm%num_casts + 1) then
            index = storm%num_casts + 1
        else
            t0 = storm%track(1,last_storm_index - 1)
            t1 = storm%track(1,last_storm_index)
            if (t0 < t .and. t <= t1) then
                index = last_storm_index
            else if ( t1 < t ) then
                found = .false.
                do index=last_storm_index+1,storm%num_casts
                    if (t < storm%track(1,index)) then
                        found = .true.
                        exit
                    endif
                enddo
                ! Assume we have gone past last forecast time
                if (.not. found) then
                    index = storm%num_casts + 1
                endif
            else ! t <= t0
                if (last_storm_index == 2) then
                    index = -1
                else
                    do index=last_storm_index-1,2,-1
                        if (storm%track(1,index-1) < t) exit
                    enddo
                endif
            endif
        endif

    end function storm_index


    ! ==========================================================================
    !  get_silva_storm_data()
    !    Interpolates in time and returns storm data.
    ! ==========================================================================
    subroutine get_silva_storm_data(t, storm, location, velocity, &
                                                max_wind_radius,    &
                                                max_wind_speed,     &
                                                central_pressure,   &
                                                rrp)

        use geoclaw_module, only: deg2rad, latlon2xy, xy2latlon

        implicit none

        ! Input
        real(kind=8), intent(in) :: t                       ! Current time
        type(silva_storm_type), intent(in) :: storm   ! Storm

        ! Output
        real(kind=8), intent(out) :: location(2), velocity(2)
        real(kind=8), intent(out) :: max_wind_radius, max_wind_speed
        real(kind=8), intent(out) :: central_pressure, rrp

        ! Local
        real(kind=8) :: fn(8), fnm(8), weight, tn, tnm, x(2)
        integer :: i

        ! Increment storm data index if needed and not at end of forecast
        i = storm_index(t,storm)
        last_storm_index = i

        ! List of possible error conditions
        if (i <= 1) then
            if (i == 0) then
                print *,"Invalid storm forecast requested for t = ",t
                print *,"Time requested is before any forecast data."
                print *,"    first time = ",storm%track(1,1)
                print *,"   ERROR = ",i
                stop
            else if (i > storm%num_casts + 2) then
                print *,"Invalid storm indexing, i > num_casts + 2..."
                print *,"This really should not happen, what have you done?"
                stop
            endif
        endif

        ! Interpolate in time for all parameters
        if (i == storm%num_casts + 1) then
            i = i - 1
            ! At last forecast, use last data for storm strength parameters and
            ! velocity, location uses last velocity and constant motion forward

            ! Convert coordinates temporarily to meters so that we can use
            ! the pre-calculated m/s velocities from before
            x = latlon2xy(storm%track(2:3,i),storm%track(2:3,i))
            x = x + (t - storm%track(1,i)) * storm%velocity(:,i)

            fn = [xy2latlon(x,storm%track(2:3,i)), &
                  storm%velocity(:,i), storm%max_wind_radius(i), &
                  storm%max_wind_speed(i), storm%central_pressure(i), &
                  storm%rrp(i)]
        else
            ! Inbetween two forecast time points (the function storm_index
            ! ensures that we are not before the first data point, i.e. i > 1)
            tn = storm%track(1,i)
            tnm = storm%track(1,i-1)
            weight = (t - tnm) / (tn - tnm)
            fn = [storm%track(2:3,i),storm%velocity(:,i), &
                  storm%max_wind_radius(i),storm%max_wind_speed(i), &
                  storm%central_pressure(i), storm%rrp(i)]
            fnm = [storm%track(2:3,i - 1),storm%velocity(:,i - 1), &
                   storm%max_wind_radius(i - 1),storm%max_wind_speed(i - 1), &
                  storm%central_pressure(i - 1), storm%rrp(i - 1)]
            fn = weight * (fn - fnm) + fnm
        endif

        ! Set output variables
        location = fn(1:2)
        velocity = fn(3:4)
        max_wind_radius = fn(5)
        max_wind_speed = fn(6)
        central_pressure = fn(7)
        rrp = fn(8)

    end subroutine get_silva_storm_data


    ! ==========================================================================
    !  set_silva_storm_fields()
    ! ==========================================================================
    subroutine set_silva_storm_fields(maux,mbc,mx,my,xlower, &
                                    ylower,dx,dy,t,aux, wind_index, &
                                    pressure_index, storm)

        use geoclaw_module, only: g => grav, rho_air, ambient_pressure
        use geoclaw_module, only: coriolis, deg2rad
        use geoclaw_module, only: spherical_distance

        use geoclaw_module, only: rad2deg

        implicit none

        ! Time of the wind field requested
        integer, intent(in) :: maux,mbc,mx,my
        real(kind=8), intent(in) :: xlower,ylower,dx,dy,t

        ! Storm description, need in out here since we may update the storm
        ! if at next time point
        type(silva_storm_type), intent(in out) :: storm

        ! Array storing wind and pressure field
        integer, intent(in) :: wind_index, pressure_index
        real(kind=8), intent(inout) :: aux(maux,1-mbc:mx+mbc,1-mbc:my+mbc)

        ! Local storage
        real(kind=8) :: x, y, r, theta, sloc(2)
        real(kind=8) :: f, mwr, mws, Pc, Pa, dp, wind, tv(2), rrp
        real(kind=8) :: mod_mws, trans_speed, ramp
        integer :: i,j
        real(kind=8) :: Fv, Nc, teta, alpha


        ! Get storm data interpolated in time
        call get_silva_storm_data(t,storm,sloc,tv,mwr,mws,Pc,rrp)

        !!! I dati dello strom sono in S.I.
        !!! ERN fa i calcoli in altre u.d.m.
        !!! converto qui e poi converno indietroa al S.I.

        ! Other quantities of interest
        Pa = ambient_pressure

        !conversions
        Pa = Pa*0.01
        ! converti Pa to mbar
        Pc = Pc * 0.01
        ! converti m to km
        mwr = mwr * 0.001
        ! converti km/h to m/s
        mws = mws / 3.6


        ! Calculate central pressure difference
        dp = Pa - Pc
        ! Limit central pressure deficit due to bad ambient pressure,
        ! really should have better ambient pressure...
        if (dp < 1.d0) dp = 1.d0

        ! Calculate silva parameters
        ! Subtract translational speed of storm from maximum wind speed
        ! to avoid distortion in the silva curve fit.  Added back later
!         trans_speed = sqrt(tv(1)**2 + tv(2)**2)
!         mod_mws = mws - trans_speed
        ! Convert wind speed (10 m) to top of atmospheric boundary layer
!         mod_mws = mod_mws / atmos_boundary_layer

        !! UR = 21.8 *sqrt(PN − P0) − 0.5* f*R [km/h]
        !! f Coriolis should be converted in rad/h

        mod_mws = 21.8*sqrt(dp)-0.5* coriolis(sloc(2))/3600. *mwr

        ! Set initial wind and pressure field, do not really need to do this
        aux(wind_index:wind_index+1,:,:) = 0.d0
        aux(pressure_index,:,:) = Pa

        ! angle of direction of cyclon
        teta = atan2( tv(1), tv(2) )


        ! Set fields
        do j=1-mbc,my+mbc
            y = ylower + (j-0.5d0) * dy     ! Degrees latitude
            f = coriolis(y)
            Nc = f/3600. * mwr / mod_mws

            do i=1-mbc,mx+mbc
                x = xlower + (i-0.5d0) * dx   ! Degrees longitude

                ! Calculate storm centric polar coordinate location of grid
                ! cell center, uses Haversine formula
                ! r in m...
                r = spherical_distance(x, y, sloc(1), sloc(2))

                ! Set pressure field and convert [mb] -- > [Pa]
                aux(pressure_index,i,j) = (Pc + dp * exp(-(mwr*1000. / r)))*100.d0

                ! Compute angles:
                ! theta: angle of point x,y respect to cyclon eye
                !! controlla gli angoli sono sbagliati per ERN
                theta = atan2((y - sloc(2)),(x - sloc(1)))
                alpha = atan2((x - sloc(1)), (y - sloc(2)))

                ! r in m...
                call calc_Fv(Nc, LOG10(r/mwr/1000.), Fv)
                ! Speed of wind at this point
                !!! wind = 0.886* (Fv * UR + 0.5* Vf_mod * np.cos(teta-beta))
                wind = 0.886* ( Fv * mod_mws + &
                     & 0.5* SQRT(tv(1)*tv(1)+tv(2)*tv(2)) * cos(alpha-teta) )

                ! Also convert from 8 minute to 1 minute
                ! sustained winds
                ! Vc = 0.002 * Vm**2 + 0.9953 * Vm
                wind = wind* (0.0012 *wind+ 1.1114 )

                ! First convert from km/h to m/s
                ! Then convert to 10-min average wind,
                ! with conversion factor G= 1.11
                ! average wind=(max 1-min sustained wind )/ G
                !!! wind = 3.6 * wind : conversion [km/h] --> [m/s]
                wind = wind*3.6/1.11

                ! Velocity components of storm (assumes perfect vortex shape)
                aux(wind_index,i,j)   = -wind * sin(theta)
                aux(wind_index+1,i,j) =  wind * cos(theta)

                !!! Apply distance ramp down(up) to fields to limit scope
                !!ramp = 0.5d0 * (1.d0 - tanh((r - rrp) / RAMP_WIDTH))
                !!aux(pressure_index,i,j) = Pa + (aux(pressure_index,i,j) - Pa) &
                !!                        * ramp
                !!aux(wind_index:wind_index+1,i,j) =                        &
                !!                        aux(wind_index:wind_index+1,i,j)  &
                !!                        * ramp
            enddo
        enddo

    end subroutine set_silva_storm_fields


    subroutine calc_Fv(Nc, log10rR, Fv)
        real(kind=8),intent(out) :: Fv
        real(kind=8),intent(in) :: log10rR, Nc
        real(kind=8) :: a,b,c,d

        if (log10rR<=0.0)then
            a  = -0.233
            b  = -12.91
            c  = -19.38
            d  = -8.311
        else
            if(Nc<=0.005)then
                a = 0.033 +  Nc*( -16.1   + 161.9* Nc )
                b = -0.43 +  Nc*(  38.9   - 316. * Nc )
                c = 0.113 +  Nc*( -28.6   + 71.1 * Nc )
                d =       +  Nc*(   1.818 + 80.6 * Nc )
            else
                a = -0.175 + Nc*(- 0.76 + Nc*(   11.7 + Nc*(- 28.1 +  17.* Nc )))
                b =  0.235 + Nc*(  2.71 + Nc*( - 67.6 + Nc*( 189.  - 155.* Nc )))
                c = -0.468 + Nc*(- 9.   + Nc*(   87.8 + Nc*(-224.  + 183*  Nc )))
                d =  0.082 + Nc*(+ 3.33 + Nc*( - 26.  + Nc*(  63.8 -  51.4*Nc )))
            end if
        end if
        Fv = log10rR*(a+log10rR*(b+log10rR*(c+d*log10rR)))
    end subroutine
end module silva_storm_module






