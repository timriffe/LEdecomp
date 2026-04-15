

library(devtools)
load_all()
library(usethis)
library(rhub)
document()

devtools::check()
check_win_devel()      # sent 15-4-2026
check_win_release()    # sent 15-4-2026
check_win_oldrelease() # sent 15-4-2026
rhub_platforms()
rhub::rhub_check(platforms = c("linux","macos","macos-arm64","windows"))

library(spelling)
spell_check()

# doesn't work until on CRAN?
#library(revdepcheck)
#revdepcheck::revdep_check(cran = FALSE)

devtools::release()
