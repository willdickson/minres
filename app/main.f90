program main
    use minres, only : minres_ez_t
    implicit none

    type(minres_ez_t) :: minres_ez

    call minres_ez % print

end program main
