.onload <- function(libname, pkgname)
{
    library.dynam("EGAnet",package=pkgname,lib.loc=libname)
}

.onAttach <- function(libname, pkgname)
{
    msg <- styletext(styletext(paste("\nEGAnet (version ", packageVersion("EGAnet"), ")", sep = ""), defaults = "underline"), defaults = "bold")
    msg <- paste(msg,'\n\nFor help getting started, see <https://r-ega.net>')
    msg <- paste(msg,"\n\nFor bugs and errors, submit an issue to <https://github.com/hfgolino/EGAnet/issues>")

    packageStartupMessage(msg)
}
