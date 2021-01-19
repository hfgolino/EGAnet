.onload <- function(libname, pkgname)
{library.dynam("EGAnet",package=pkgname,lib.loc=libname)}

.onAttach <- function(libname, pkgname)
{
    temp <- packageDescription("EGAnet")
    msg <- paste("Package: ",temp$Package,": ",temp$Title,"\n",
               "Version: ",temp$Version,"\n",
               "Updated on: ",
               temp$Date,"\n", sep="")
    msg <- paste(msg,"Maintainer: Hudson Golino, University of Virginia\n",sep="")
    msg <- paste(msg, "Authors: Hudson Golino, University of Virginia & \n Alexander Christensen, University of Pennsylvania\n",sep="")
    msg <- paste(msg,"Contributors: Robert Moulder, University of Virginia,\n",
                 "Luis Eduardo Garrido, Pontificia Universidad Catolica Madre y Maestra\n",sep="")
    msg <- paste(msg,'For citation information, type citation("EGAnet")\n')
    msg <- paste(msg,'For vignettes, type browseVignettes("EGAnet")\n')
    msg <- paste(msg,"For bugs and errors, submit an issue to <https://github.com/hfgolino/EGAnet/issues>")
    packageStartupMessage(msg)
}
