#' @name accident2014
#' @title Sample of car accident location in the UK during year 2014.
#' @description Longitude and latitude of 500 car accident during year 2014 (source: data.gov.uk).
#' @docType data
#' @usage accident2014
#' @format The dataset has 500 instances described by 2 variables (coordinates).
#' @source \url{https://data.gov.uk/}
NULL

#' @name alcohol
#' @title Alcohol dataset
#' @description This dataset has been extracted from the WHO database and depict the alcool habits in the 27 european contries (in 2010).
#' @docType data
#' @usage alcohol
#' @format The dataset has 27 instances described by 4 variables.
#' The variables are the average amount of alcool of different types per year par inhabitent.
#' @source \url{https://www.who.int/}
NULL

#' @name autompg
#' @title Auto MPG dataset
#' @description This dataset was taken from the StatLib library which is maintained at Carnegie Mellon University.
#' The dataset was used in the 1983 American Statistical Association Exposition.
#' @docType data
#' @usage autompg
#' @format The dataset has 392 instances described by 8 variables.
#' The seven first variables are numeric variables. The last variable is qualitative (car origin).
#' @source \url{https://archive.ics.uci.edu/ml/datasets/auto+mpg}
NULL

#' @name beetles
#' @title Flea beetles dataset
#' @description Data were collected on the genus of flea beetle \emph{Chaetocnema}, which contains three species:
#' \emph{concinna}, \emph{heikertingeri}, and \emph{heptapotamica}.
#' Measurements were made on the width and angle of the aedeagus of each beetle.
#' The goal of the original study was to form a classification rule to distinguish the three species.
#' @docType data
#' @usage beetles
#' @format The dataset has 74 instances described by 3 variables.
#' The variables are as follows:
#' \describe{
#' \item{\code{Width}}{The maximal width of aedeagus in the forpart (in microns).}
#' \item{\code{Angle}}{The front angle of the aedeagus (1 unit = 7.5 degrees).}
#' \item{\code{Shot.put}}{Species of flea beetle from the genus \emph{Chaetocnema}.}
#' }
#' @source Lubischew, A.A. (1962) On the use of discriminant functions in taxonomy. Biometrics, 18, 455-477.
NULL

#' @name britpop
#' @title Population and location of 18 major british cities.
#' @description Longitude and latitude and population of 18 major cities in the Great Britain.
#' @docType data
#' @usage britpop
#' @format The dataset has 18 instances described by 3 variables.
NULL

#' @name cookies
#' @aliases cookies.desc.train cookies.desc.test cookies.y.train cookies.y.test
#' @title Cookies dataset
#' @description This data set contains measurements from quantitative NIR spectroscopy.
#' The example studied arises from an experiment done to test the feasibility of NIR spectroscopy to measure the composition of biscuit dough pieces (formed but unbaked biscuits).
#' Two similar sample sets were made up, with the standard recipe varied to provide a large range for each of the four constituents under investigation: fat, sucrose, dry flour, and water.
#' The calculated percentages of these four ingredients represent the 4 responses.
#' There are 40 samples in the calibration or training set (with sample 23 being an outlier).
#' There are a further 32 samples in the separate prediction or validation set (with example 21 considered as an outlier).
#' An NIR reflectance spectrum is available for each dough piece.
#' The spectral data consist of 700 points measured from 1100 to 2498 nanometers (nm) in steps of 2 nm.
#' @docType data
#' @usage cookies
#' cookies.desc.train
#' cookies.desc.test
#' cookies.y.train
#' cookies.y.test
#' @format The cookies.desc.* datasets contains the 700 columns that correspond to the NIR reflectance spectrum.
#' The cookies.y.* datasets contains four columns that correspond to the four constituents fat, sucrose, dry flour, and water.
#' The cookies.*.train contains 40 rows that correspond to the calibration data.
#' The cookies.*.test contains 32 rows that correspond to the prediction data.
#' @source P. J. Brown and T. Fearn and M. Vannucci (2001) "Bayesian wavelet regression on curves with applications to a spectroscopic calibration problem", Journal of the American Statistical Association, 96(454), pp. 398-408.
#' @seealso \code{\link[fds]{labp}}, \code{\link[fds]{labc}}, \code{\link[fds]{nirp}}, \code{\link[fds]{nirc}}
NULL

#' @name credit
#' @title Credit dataset
#' @description This is a fake dataset simulating a bank database about loan clients.
#' @docType data
#' @usage credit
#' @format The dataset has 66 instances described by 11 qualitative variables.
NULL

#' Parabol dataset
#'
#' Generate a random dataset shaped like a parabol and a gaussian distribution
#' @name data.parabol
#' @param n Number of observations in each class.
#' @param xlim Minimum and maximum on the x axis.
#' @param center Coordinates of the center of the gaussian distribution.
#' @param coeff Coefficient of the parabol.
#' @param sigma Variance in each class.
#' @param levels Name of each class.
#' @param graph A logical indicating whether or not a graphic should be plotted.
#' @param seed A specified seed for random number generation.
#' @return A randomly generated dataset.
#' @export
#' @seealso \code{\link{data.target1}}, \code{\link{data.target2}}, \code{\link{data.twomoons}}
#' @examples
#' data.parabol ()
data.parabol <-
  function (n = c (500, 100), xlim = c (-3, 3), center = c (0, 4), coeff = 0.5, sigma = c (0.5, 0.5), levels = NULL, graph = TRUE, seed = NULL)
  {
    set.seed (seed)
    if (is.null (levels))
      levels = paste ("Class", 1:2)
    d = stats::runif (n [1], xlim [1], xlim [2])
    d = cbind (d, coeff * d^2 + stats::rnorm (n [1], 0, sigma [1]))
    d = rbind (d, cbind (stats::rnorm (n [2], center [1], sigma [2]), stats::rnorm (n [2], center [2], sigma [2])))
    d = cbind.data.frame (d, factor (c (rep (1, n [1]), rep (2, n [2])), labels = levels))
    colnames (d) = c ("X", "Y", "Class")
    if (graph)
      plotdata (d [, -3], d [, 3])
    return (d)
  }

#' Target1 dataset
#'
#' Generate a random dataset shaped like a target.
#' @name data.target1
#' @param r Radius of each class.
#' @param n Number of observations in each class.
#' @param sigma Variance in each class.
#' @param levels Name of each class.
#' @param graph A logical indicating whether or not a graphic should be plotted.
#' @param seed A specified seed for random number generation.
#' @return A randomly generated dataset.
#' @export
#' @seealso \code{\link{data.parabol}}, \code{\link{data.target2}}, \code{\link{data.twomoons}}
#' @examples
#' data.target1 ()
data.target1 <-
  function (r = 1:3, n = 200, sigma = .1, levels = NULL, graph = TRUE, seed = NULL)
  {
    set.seed (seed)
    if (length (n) == 1)
      n = rep (n, length (r))
    if (is.null (levels))
      levels = paste ("Class", 1:length (r))
    alpha = stats::runif (sum (n), 0, 2 * pi)
    k = as.vector (sapply (1:length (r), function (index) rep (r [index], n [index])))
    r = stats::rnorm (sum (n), 0, .1) + k
    x = r * cos (alpha)
    y = r * sin (alpha)
    d = cbind.data.frame (x, y, factor (k, labels = levels))
    colnames (d) = c ("X", "Y", "Class")
    if (graph)
      plotdata (d [, -3], d [, 3])
    return (d)
  }


#' Target2 dataset
#'
#' Generate a random dataset shaped like a target.
#' @name data.target2
#' @param minr Minimum radius of each class.
#' @param maxr Maximum radius of each class.
#' @param initn Number of observations at the beginning of the generation process.
#' @param levels Name of each class.
#' @param graph A logical indicating whether or not a graphic should be plotted.
#' @param seed A specified seed for random number generation.
#' @return A randomly generated dataset.
#' @export
#' @seealso \code{\link{data.parabol}}, \code{\link{data.target1}}, \code{\link{data.twomoons}}
#' @examples
#' data.target2 ()
data.target2 <-
  function (minr = c (0, 2), maxr = minr + 1, initn = 1000, levels = NULL, graph = TRUE, seed = NULL)
  {
    set.seed (seed)
    limits = max (maxr)
    if (is.null (levels))
      levels = paste ("Class", 1:2)
    d = matrix (stats::runif (2 * initn, -limits, limits), ncol = 2)
    d = cbind (d, t (flexclust::dist2 (matrix (c (0, 0), ncol = 2), d)))
    d = d [d [, 3] <= limits, ]
    d = cbind.data.frame (d, factor (sapply (d [, 3], function (dis) min (which (dis <= maxr))), labels = levels))
    d = d [apply (sapply (1:length (minr), function (index) return ((d [, 3] >= minr [index]) & (d [, 3] <= maxr [index]))), 1, any), -3]
    colnames (d) = c ("X", "Y", "Class")
    if (graph)
      plotdata (d [, -3], d [, 3])
    return (d)
  }

#' Two moons dataset
#'
#' Generate a random dataset shaped like two moons.
#' @name data.twomoons
#' @param r Radius of each class.
#' @param n Number of observations in each class.
#' @param sigma Variance in each class.
#' @param levels Name of each class.
#' @param graph A logical indicating whether or not a graphic should be plotted.
#' @param seed A specified seed for random number generation.
#' @return A randomly generated dataset.
#' @export
#' @seealso \code{\link{data.parabol}}, \code{\link{data.target1}}, \code{\link{data.target2}}
#' @examples
#' data.twomoons ()
data.twomoons <-
  function (r = 1, n = 200, sigma = .1, levels = NULL, graph = TRUE, seed = NULL)
  {
    set.seed (seed)
    if (length (n) == 1)
      n = rep (n, 2)
    if (is.null (levels))
      levels = paste ("Class", 1:2)
    alpha1 = stats::runif (n [1], 0, pi)
    alpha2 = stats::runif (n [2], 0, pi)
    k = c (rep (1, n [1]), rep (2, n [2]))
    sx = r / 2
    sy = r / 6
    x = c (r * cos (alpha1) + sx, r * cos (alpha2) - sx)
    y = c (r * sin (alpha1) - sy, r * -sin (alpha2) + sy)
    d = cbind.data.frame (x, y, factor (k, labels = levels))
    noise = matrix (stats::rnorm (2 * sum (n), 0, sigma), ncol = 2)
    d  [, 1:2] = d  [, 1:2] + noise
    colnames (d) = c ("X", "Y", "Class")
    if (graph)
      plotdata (d [, -3], d [, 3])
    return (d)
  }

#' @name data1
#' @title "data1" dataset
#' @description Synthetic dataset.
#' @docType data
#' @usage data1
#' @format 240 observations described by 4 variables and grouped into 16 classes.
#' @author Alexandre Blansché \email{alexandre.blansche@univ-lorraine.fr}
NULL

#' @name data2
#' @title "data2" dataset
#' @description Synthetic dataset.
#' @docType data
#' @usage data2
#' @format 500 observations described by 10 variables and grouped into 3 classes.
#' @author Alexandre Blansché \email{alexandre.blansche@univ-lorraine.fr}
NULL

#' @name data3
#' @title "data3" dataset
#' @description Synthetic dataset.
#' @docType data
#' @usage data3
#' @format 300 observations described by 3 variables and grouped into 3 classes.
#' @author Alexandre Blansché \email{alexandre.blansche@univ-lorraine.fr}
NULL

#' @name decathlon
#' @title Decathlon dataset
#' @description The dataset contains results from two athletics competitions.
#' The 2004 Olympic Games in Athens and the 2004 Decastar.
#' @docType data
#' @usage decathlon
#' @format The dataset has 41 instances described by 13 variables.
#' The variables are as follows:
#' \describe{
#' \item{\code{100m}}{In seconds.}
#' \item{\code{Long.jump}}{In meters.}
#' \item{\code{Shot.put}}{In meters.}
#' \item{\code{High.jump}}{In meters.}
#' \item{\code{400m}}{In seconds.}
#' \item{\code{110m.h}}{In seconds.}
#' \item{\code{Discus.throw}}{In meters.}
#' \item{\code{Pole.vault}}{In meters.}
#' \item{\code{Javelin.throw}}{In meters.}
#' \item{\code{1500m}}{In seconds.}
#' \item{\code{Rank}}{The rank at the competition.}
#' \item{\code{Points}}{The number of points obtained by the athlete.}
#' \item{\code{Competition}}{\code{Olympics} or \code{Decastar}.}
#' }
#' @source \url{https://husson.github.io/data.html}
NULL

#' @name eucalyptus
#' @title Eucalyptus dataset
#' @description Measuring the height of a tree is not an easy task. Is it possible to estimate the height as a function of the circumference of the trunk?
#' @docType data
#' @usage eucalyptus
#' @format The dataset has 1429 instances (eucalyptus trees) with 2 measurements: the height and the circumference.
#' @source \url{http://www.cmap.polytechnique.fr/~lepennec/enseignement/MAP553/Lab2_Linear.html}
NULL

#' @name ionosphere
#' @title Ionosphere dataset
#' @description This is a dataset from the UCI repository.
#' This radar data was collected by a system in Goose Bay, Labrador. This system consists of a phased array of 16 high-frequency antennas with a total transmitted power on the order of 6.4 kilowatts. See the paper for more details. The targets were free electrons in the ionosphere. "Good" radar returns are those showing evidence of some type of structure in the ionosphere. "Bad" returns are those that do not; their signals pass through the ionosphere.
#' Received signals were processed using an autocorrelation function whose arguments are the time of a pulse and the pulse number. There were 17 pulse numbers for the Goose Bay system. Instances in this databse are described by 2 attributes per pulse number, corresponding to the complex values returned by the function resulting from the complex electromagnetic signal.
#' One attribute with constant value has been removed.
#' @docType data
#' @usage ionosphere
#' @format The dataset has 351 instances described by 34. The last variable is the class.
#' @source \url{https://archive.ics.uci.edu/ml/datasets/ionosphere}
NULL

#' @name linsep
#' @title Linsep dataset
#' @description Synthetic dataset.
#' @docType data
#' @usage linsep
#' @format Class \code{A} contains 50 observations and class \code{B} contains 500 observations.
#' There are two numeric variables: \code{X} and \code{Y}.
#' @author Alexandre Blansché \email{alexandre.blansche@univ-lorraine.fr}
NULL

#' @name movies
#' @title Movies dataset
#' @description Extract from the movie lens dataset. Missing values have been imputed.
#' @docType data
#' @usage movies
#' @format A set of 49 movies, rated by 55 users.
#' @source \url{https://grouplens.org/datasets/movielens/}
NULL

#' @name ozone
#' @title Ozone dataset
#' @description This dataset constains measurements on ozone level.
#' @docType data
#' @usage ozone
#' @format Each instance is described by the maximum level of ozone measured during the day.
#' Temperature, clouds, and wind are also recorded.
#' @source \url{https://r-stat-sc-donnees.github.io/ozone.txt}
NULL

#' @name reg1
#' @aliases reg1.train reg1.test
#' @title reg1 dataset
#' @description Artificial dataset for simple regression tasks.
#' @docType data
#' @usage reg1
#' reg1.train
#' reg1.test
#' @format 50 instances and 3 variables. \code{X}, a numeric, \code{K}, a factor, and \code{Y}, a numeric (the target variable).
#' @author Alexandre Blansché \email{alexandre.blansche@univ-lorraine.fr}
NULL

#' @name reg2
#' @aliases reg2.train reg2.test
#' @title reg2 dataset
#' @description Artificial dataset for simple regression tasks.
#' @docType data
#' @usage reg2
#' reg2.train
#' reg2.test
#' @format 50 instances and 2 variables. \code{X} and \code{Y} (the target variable) are both numeric variables.
#' @author Alexandre Blansché \email{alexandre.blansche@univ-lorraine.fr}
NULL

#' @name snore
#' @title Snore dataset
#' @description This dataset has been used in a study on snoring in Angers hospital.
#' @docType data
#' @usage snore
#' @format The dataset has 100 instances described by 7 variables.
#' The variables are as follows:
#' \describe{
#' \item{\code{Age}}{In years.}
#' \item{\code{Weights}}{In kg.}
#' \item{\code{Height}}{In cm.}
#' \item{\code{Alcool}}{Number of glass of alcool per day.}
#' \item{\code{Sex}}{M for male or F for female.}
#' \item{\code{Snore}}{Snoring diagnosis (Y or N).}
#' \item{\code{Tobacco}}{Y or N.}
#' }
#' @source \url{http://forge.info.univ-angers.fr/~gh/Datasets/datasets.htm}
NULL

#' @name spine
#' @aliases spine.train spine.test
#' @title Spine dataset
#' @description The data have been organized in two different but related classification tasks.
#' The first task consists in classifying patients as belonging to one out of three categories: Normal, Disk Hernia or Spondylolisthesis.
#' For the second task, the categories Disk Hernia and Spondylolisthesis were merged into a single category labelled as 'abnormal'.
#' Thus, the second task consists in classifying patients as belonging to one out of two categories: Normal or Abnormal.
#' @docType data
#' @usage spine
#' spine.train
#' spine.test
#' @format The dataset has 310 instances described by 8 variables.
#' Variables V1 to V6 are biomechanical attributes derived from the shape and orientation of the pelvis and lumbar spine.
#' The variable Classif2 is the classification into two classes \code{AB} and \code{NO}.
#' The variable Classif3 is the classification into 3 classes \code{DH}, \code{SL} and \code{NO}.
#' \code{spine.train} contains 217 instances and \code{spine.test} contains 93.
#' @source \url{http://archive.ics.uci.edu/ml/datasets/vertebral+column}
NULL

#' @name temperature
#' @title Temperature dataset
#' @description The data contains temperature measurement and geographic coordinates of 35 european cities.
#' @docType data
#' @usage temperature
#' @format The dataset has 35 instances described by 17 variables.
#' Average temperature of the 12 month. Mean and amplitude of the temperature. Latitude and longitude of the city. Localisation in Europe.
NULL

#' @name titanic
#' @title Titanic dataset
#' @description This dataset from the British Board of Trade depict the fate of the passengers and crew during the RMS Titanic disaster.
#' @docType data
#' @usage titanic
#' @format The dataset has 2201 instances described by 4 variables.
#' The variables are as follows:
#' \describe{
#' \item{\code{Category}}{1st, 2nd, 3rd Class or Crew.}
#' \item{\code{Age}}{Adult or Child.}
#' \item{\code{Sex}}{Female or Male.}
#' \item{\code{Fate}}{Casualty or Survivor.}
#' }
#' @source British Board of Trade (1990), Report on the Loss of the ‘Titanic’ (S.S.). British Board of Trade Inquiry Report (reprint). Gloucester, UK: Allan Sutton Publishing.
#' @seealso \code{\link[datasets]{Titanic}}
NULL

#' @name birth
#' @title Birth dataset
#' @description Tutorial data set (vector).
#' @docType data
#' @usage birth
#' @format The dataset is a names vector of nine values (birth years).
NULL

#' @name universite
#' @title University dataset
#' @description The dataset presents a french university demographics.
#' @docType data
#' @usage universite
#' @format The dataset has 10 instances (university departments) described by 12 variables.
#' The fist six variables are the number of female and male student
#' studying for bachelor degree (Licence), master degree (Master) and doctorate (Doctorat).
#' The six last variables are obtained by combining the first ones.
#' @source \url{https://husson.github.io/data.html}
NULL

#' @name wine
#' @title Wine dataset
#' @description These data are the results of a chemical analysis of wines grown in the same region in Italy but derived from three different cultivars.
#' The analysis determined the quantities of 13 constituents found in each of the three types of wines.
#' @docType data
#' @usage wine
#' @format There are 178 observations and 14 variables.
#' The first variable is the class label (\code{1}, \code{2}, \code{3}).
#' @source \url{https://archive.ics.uci.edu/ml/datasets/wine}
NULL

#' @name wheat
#' @title Wheat dataset
#' @description The data contains kernels belonging to three different varieties of wheat: Kama, Rosa and Canadian, 70 elements each, randomly selected.
#' High quality visualization of the internal kernel structure was detected using a soft X-ray technique. The images were recorded on 13x18 cm X-ray KODAK plates.
#' Source : Institute of Agrophysics of the Polish Academy of Sciences in Lublin.
#' @docType data
#' @usage wheat
#' @format The dataset has 210 instances described by 8 variables:
#' area, perimeter, compactness, length, width, asymmetry coefficient, groove length and variery.
#' @source \url{https://archive.ics.uci.edu/ml/datasets/seeds}
NULL

#' @name zoo
#' @title Zoo dataset
#' @description Animal description based on various features.
#' @docType data
#' @usage zoo
#' @format The dataset has 101 instances described by 17 qualitative variables.
#' @source \url{https://archive.ics.uci.edu/ml/datasets/zoo}
NULL
