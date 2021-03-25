.onload <- function(libname, pkgname)
{library.dynam("EGAnet",package=pkgname,lib.loc=libname)}

.onAttach <- function(libname, pkgname)
{
    msg <- styletext(styletext(paste("\nEGAnet (version ", packageVersion("EGAnet"), ")", sep = ""), defaults = "underline"), defaults = "bold")

    msg <- paste(msg,'\nFor help getting started, type browseVignettes("EGAnet")\n')
    msg <- paste(msg,"For bugs and errors, submit an issue to <https://github.com/hfgolino/EGAnet/issues>")

    msg <- paste(msg, paste(styletext("\n\nNEW", defaults = "bold"),
                            ": EGAnet will write your Methods section for you. Type ?methods.section for more details",
                            sep = ""))
    
    msg <- paste(msg, paste(styletext("\n\nMAJOR CHANGE", defaults = "bold"),
                            ": Unidimensionality adjustment in EGA functions has changed to the leading eigenvalue adjustment (see Christensen, Garrido, & Golino, 2021)",
                            sep = ""))
    packageStartupMessage(msg)
}
