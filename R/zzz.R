.onload <- function(libname, pkgname)
{
    library.dynam("EGAnet",package=pkgname,lib.loc=libname)
}

.onAttach <- function(libname, pkgname)
{
    msg <- styletext(styletext(paste("\nEGAnet (version ", packageVersion("EGAnet"), ")", sep = ""), defaults = "underline"), defaults = "bold")
    msg <- paste(msg,'\nCheck out what\'s changed: <https://tinyurl.com/EGAnet-changes>\n')
    msg <- paste(msg,"\nFor bugs and errors, submit an issue to <https://github.com/hfgolino/EGAnet/issues>")

    packageStartupMessage(msg)
}
