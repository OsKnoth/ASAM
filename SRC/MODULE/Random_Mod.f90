module numz
  integer, parameter:: b8 = selected_real_kind(14)
  real(b8), parameter :: pi = 3.141592653589793239_b8
  integer gene_size
end module 
MODULE Random_Mod
! module contains three functions
! ran1 returns a uniform random number between 0-1
! spread returns random number between min - max
! normal returns a normal distribution

contains
    function ran1()  !returns random number between 0 - 1
        use numz
        implicit none
        real(b8) ran1,x
        call init_random_seed
        call random_number(x) ! built in fortran 90 random number function
        ran1=x
    end function ran1

    function ran_spread(min,max)  !returns random number between min - max
        use numz
        implicit none
        real(b8) ran_spread
        real(b8) min,max
        ran_spread=(max - min) * ran1() + min
    end function ran_spread

    function normal(mean,sigma) !returns a normal distribution
        use numz
        implicit none
        real(b8) normal,tmp
        real(b8) mean,sigma
        integer flag
        real(b8) fac,gsave,rsq,r1,r2
        save flag,gsave
        data flag /0/
        if (flag.eq.0) then
        rsq=2.0_b8
            do while(rsq.ge.1.0_b8.or.rsq.eq.0.0_b8) ! new from for do
                r1=2.0_b8*ran1()-1.0_b8
                r2=2.0_b8*ran1()-1.0_b8
                rsq=r1*r1+r2*r2
            enddo
            fac=sqrt(-2.0_b8*log(rsq)/rsq)
            gsave=r1*fac
            tmp=r2*fac
            flag=1
        else
            tmp=gsave
            flag=0
        endif
        normal=tmp*sigma+mean
        return
    end function normal

    SUBROUTINE init_random_seed
      INTEGER :: i, n, clock
      INTEGER, DIMENSION(:), ALLOCATABLE :: seed
      INTEGER, SAVE :: inc=0

      CALL RANDOM_SEED(size = n)
      ALLOCATE(seed(n))

      CALL SYSTEM_CLOCK(COUNT=clock)

      seed = clock + 37 * (/ (i - 1, i = 1, n) /) + inc
      CALL RANDOM_SEED(PUT = seed)

      DEALLOCATE(seed)
      
      inc = inc + 1
    END SUBROUTINE

end module Random_Mod
