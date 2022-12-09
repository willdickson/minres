module minres

    use minres_mod,     only : minres_solver
    use minres_ez_mod,  only : minres_ez_t

    implicit none
    private

    public :: minres_solver
    public :: minres_ez_t

end module minres
