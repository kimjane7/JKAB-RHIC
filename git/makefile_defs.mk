#edit these to your heart's content
MADAI_HOME = /home/kimjane7/janegit
#prefix of .../rhic directory (root of madai installation)

MADAI_COMPILERDEFS =
MADAI_GSLPATH = /usr/local
#where gnu scientific library is installed

MADAI_X11PATH = /usr/X11R6
#root of X11 installation

MADAI_EIGEN_HOME = /usr/local

#MADAI_CPP = /usr/bin/clang++
MADAI_CPP = /usr/bin/g++
#compiler

MADAI_BASELIBS = -lcoralutils -lgsl -lgslcblas -lrhicutils

MADAI_INSTALLDIR = ~/local
#location of where you want things installed

#MADAI_CFLAGS = -O
#MADAI_CFLAGS = -fast
#MADAI_CFLAGS = -O2
MADAI_CFLAGS = -O2 -std=c++11
#compiler optimization flags, usually -O2 for linux, -fast for OSX with g++
