
#LIB=/home/xiaywang/pulp/PULP-DSP/lib/libmyplpdsp.a

IDIR=../../../../../include

PULP_APP = test

PULP_APP_CL_SRCS = cluster.c
PULP_APP_FC_SRCS = test.c

PULP_LDFLAGS += -lplpdsp
#PULP_LDFLAGS += $(LIB)
PULP_CFLAGS += -I$(IDIR) -O3 -g

include $(PULP_SDK_HOME)/install/rules/pulp_rt.mk
