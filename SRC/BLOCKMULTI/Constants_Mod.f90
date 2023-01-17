MODULE Constants_Mod

USE Kind_Mod


REAL(Realkind), PARAMETER :: cpd = 1004.0, &
                             cvd = 717.0, &
                             Rd = cpd - cvd, &
                             kappa = Rd / cpd, &
                             p0 = 1.0e5

END MODULE Constants_Mod
