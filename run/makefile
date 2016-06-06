include /home/kimjane7/janegit/makefile_defs.mk

RHIC_HOME = ${MADAI_HOME}/rhic

PROGDATA=${MADAI_INSTALLDIR}/progdata

INCLUDE = -I ${MADAI_INSTALLDIR}/include/madai -I ${MADAI_GSLPATH}/include

LIBRARY = -L ${MADAI_GSLPATH}/lib -L ${MADAI_INSTALLDIR}/lib/madai ${MADAI_BASELIBS}

corallibs :
	make -C ${RHIC_HOME}/coral install
	
b3dlibs :
	make -C ${RHIC_HOME}/b3d2 install

b3d : progdata_link b3dmain.cc b3dlibs corallibs
	${MADAI_CPP} -o b3d b3dmain.cc ${MADAI_CFLAGS} ${INCLUDE} ${LIBRARY} -lb3d -lcoralutils

progdata_link :
	ln -s -f ${PROGDATA} .

joshconvert : joshmain.cc b3dlibs
	${MADAI_CPP} -o joshconvert joshmain.cc ${MADAI_CFLAGS} ${INCLUDE} ${LIBRARY} -lb3d