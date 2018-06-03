library(sos)
library(htmlTable)
library(stringr)
library(dplyr)

# Get auc functions
auc.search <- findFn("auc") 
auc.functions <- auc.search %>%
	filter(Function == "auc", Package != "pROC") %>%
	select(Package, Function, Description, Link)
rownames(auc.functions) <- auc.functions$Package

# Get roc functions
roc.search <- findFn("roc") 
roc.functions <- roc.search %>%
	filter(Function == "roc", Package != "pROC") %>%
	select(Package, Function, Description, Link)
rownames(roc.functions) <- roc.functions$Package


# Install missing packages
missing.packages <- auc.functions$Package[ ! auc.functions$Package %in% installed.packages()[,"Package"]]
install.packages(missing.packages)
missing.packages <- roc.functions$Package[ ! roc.functions$Package %in% installed.packages()[,"Package"]]
install.packages(missing.packages)

# Filter packages that are still missing
available.packages.with.auc <- auc.functions[auc.functions$Package %in% installed.packages()[,"Package"],]
available.packages.with.roc <- roc.functions[roc.functions$Package %in% installed.packages()[,"Package"],]

#' Check if a function within a package is a generic function
#' @param pkg package name as a character string
#' @param fun function name as a character string
#' @return \code{TRUE} if the function is generic, \code{FALSE} otherwise. 
#' If the package doesn't contain a function named `fun`, \code{NA} is returned.
is.function.in.package.generic <- function(pkg, fun) {
	old.search.pos <- search()[2]
	on.exit({
		while (attr(parent.env(.GlobalEnv), "name") != old.search.pos) {
			detach()
		}
	})
	suppressPackageStartupMessages(library(pkg, character.only = TRUE))
	# Does the package actually have a roc function
	t <- try(get(fun), silent=TRUE)
	if (methods::is(t, "try-error")) {
		warning(sprintf("Package %s doesn't seem to contain function %s", pkg, fun))
		return(NA)
	}
	if (utils::isS3stdGeneric(fun)) {
		return(TRUE)
	}
	if (methods::isGeneric(fun)) {
		return(TRUE)
	}
	return(FALSE)
}

# Test which packages have generic functions
generics.auc <- sapply(available.packages.with.auc$Package, is.function.in.package.generic, fun="auc")
generics.roc <- sapply(available.packages.with.roc$Package, is.function.in.package.generic, fun="roc")

# Prepare table
available.packages.with.auc$Generic <- c("TRUE"="Generic", "FALSE"="Not Generic")[as.character(generics.auc)]
available.packages.with.auc$auc <- sprintf('<a href="%s">%s</a>', available.packages.with.auc$Link, available.packages.with.auc$Generic)

available.packages.with.roc$Generic <- c("TRUE"="Generic", "FALSE"="Not Generic")[as.character(generics.roc)]
available.packages.with.roc$roc <- sprintf('<a href="%s">%s</a>', available.packages.with.roc$Link, available.packages.with.roc$Generic)

# Final table
table <- data.frame(
	Package = union(available.packages.with.roc$Package, 
					available.packages.with.auc$Package))
rownames(table) <- table$Package
table[available.packages.with.roc$Package, "roc"] <- available.packages.with.roc$roc
table[available.packages.with.auc$Package, "auc"] <- available.packages.with.auc$auc

# Format as HTML table
htmlTable(table[order(table$Package), c("Package", "roc", "auc")], escape.html = FALSE)
