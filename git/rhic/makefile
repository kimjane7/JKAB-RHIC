include ../makefile_defs.mk

CORAL_PATH = rhic/coral
B3D_PATH = rhic/b3d2
RHICUTILS_PATH = rhic/rhicutils

all: coral b3d rhicutils

install: coral_install rhicutils_install b3d_install

clean: coral_clean b3d_clean rhicutils_clean

uninstall: coral_uninstall b3d_uninstall rhicutils_uninstall

coral:
	${MAKE} -C ${MADAI_HOME}/${CORAL_PATH} -f ${MADAI_HOME}/${CORAL_PATH}/makefile coral;

coral_install:
	${MAKE} -C ${MADAI_HOME}/${CORAL_PATH} -f ${MADAI_HOME}/${CORAL_PATH}/makefile install;

coral_clean:
	${MAKE} -C ${MADAI_HOME}/${CORAL_PATH} -f ${MADAI_HOME}/${CORAL_PATH}/makefile clean;

coral_uninstall:
	${MAKE} -C ${MADAI_HOME}/${CORAL_PATH} -f ${MADAI_HOME}/${CORAL_PATH}/makefile uninstall;

b3d: coral coral_install
	${MAKE} -C ${MADAI_HOME}/${B3D_PATH} -f ${MADAI_HOME}/${B3D_PATH}/makefile all;

b3d_install:
	${MAKE} -C ${MADAI_HOME}/${B3D_PATH} -f ${MADAI_HOME}/${B3D_PATH}/makefile install;

b3d_clean:
	${MAKE} -C ${MADAI_HOME}/${B3D_PATH} -f ${MADAI_HOME}/${B3D_PATH}/makefile clean;

b3d_uninstall:
	${MAKE} -C ${MADAI_HOME}/${B3D_PATH} -f ${MADAI_HOME}/${B3D_PATH}/makefile uninstall;

rhicutils: coral coral_install
	${MAKE} -C ${MADAI_HOME}/${RHICUTILS_PATH} -f ${MADAI_HOME}/${RHICUTILS_PATH}/makefile all;

rhicutils_install:
	${MAKE} -C ${MADAI_HOME}/${RHICUTILS_PATH} -f ${MADAI_HOME}/${RHICUTILS_PATH}/makefile install;

rhicutils_clean:
	${MAKE} -C ${MADAI_HOME}/${RHICUTILS_PATH} -f ${MADAI_HOME}/${RHICUTILS_PATH}/makefile clean;

rhicutils_uninstall:
	${MAKE} -C ${MADAI_HOME}/${RHICUTILS_PATH} -f ${MADAI_HOME}/${RHICUTILS_PATH}/makefile uninstall;

