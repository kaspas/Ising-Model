# ####################################################################
# 
#			   C/C++ Makefile
# 
# Adapted from
#  D.Floros' template and
#  http://www.cs.swarthmore.edu/~newhall/unixhelp/howto_makefiles.html
#  
# ####################################################################
# 
# 'make'        build executable file 'main'
# 'make lib'	build the libraries .a
# 'make clean'  removes all .o and executable files
#

# define the shell to bash
SHELL := /bin/bash

# define the C/C++ compiler to use,default here is clang
CC = gcc-7 -std=gnu11 -pg
NVCC =nvcc 
# define compile-time flags
CFLAGS = -O3
CUDACFLAGS= -O3 --use_fast_math -lineinfo -arch=sm_70 -Xptxas=-v
# define any directories containing header files
INCLUDES = -I./inc

# define library paths in addition to /usr/lib
#   if I wanted to include libraries not in /usr/lib specify
#   their path using -Lpath, something like:
LDFLAGS =

# define any libraries to link into executable
LIBS =  

# define the source file for the library
SRC = ising

# define directories for source files and libraries
SRCDIR = src
LIBDIR = lib

# define vpath to look for missing source files
VPATH = src

# define the different possible executables
TYPES = sequential 

EXTRAS = V1 V2 V3

# define the executable file  name
MAIN = main
TEST = test
S_TEST = test_rand


# call everytime
.PRECIOUS: %.a

all: $(addprefix $(MAIN)_, $(TYPES)) $(addsuffix $(MAIN), $(EXTRAS))

lib: $(addsuffix .a, $(addprefix $(LIBDIR)/$(SRC)_, $(TYPES))) $(addsuffix .a, $(addprefix $(LIBDIR)/$(SRC)_, $(EXTRAS)))

$(MAIN)_%: $(MAIN).c  $(LIBDIR)/$(SRC)_%.a 
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $(word 1,$^) $(word 2,$^) $(LDFLAGS) $(LIBS)

%$(Extras)$(MAIN): $(MAIN).c $(LIBDIR)/$(SRC)_%.a 
	cp $(SRCDIR)/$(MAIN).c $(SRCDIR)/$(MAIN).cu
	$(NVCC) $(CUDACFLAGS) $(INCLUDES) -o $@ $(word 1,$^)u $(word 2,$^) $(LDFLAGS) $(LIBS)
	rm $(SRCDIR)/$(MAIN).cu
# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .cpp file) and $@: the name of the target of the rule (a .o file)

run: all
	./main_sequential $(n) $(k)  
	./V1main $(n) $(k)
	./V2main $(n) $(k)
	./V3main $(n) $(k)


$(TEST)_%: $(TEST).c  $(LIBDIR)/$(SRC)_%.a 
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $(word 1,$^) $(word 2,$^) $(LDFLAGS) $(LIBS)

%$(Extras)$(TEST): $(TEST).c $(LIBDIR)/$(SRC)_%.a 
	cp $(SRCDIR)/$(TEST).c $(SRCDIR)/$(TEST).cu
	$(NVCC) $(CUDACFLAGS) $(INCLUDES) -o $@ $(word 1,$^)u $(word 2,$^) $(LDFLAGS) $(LIBS)
	rm $(SRCDIR)/$(TEST).cu



%$(Extras)$(S_TEST): $(S_TEST).c $(LIBDIR)/$(SRC)_%.a 
	cp $(SRCDIR)/$(S_TEST).c $(SRCDIR)/$(S_TEST).cu
	$(NVCC) $(CUDACFLAGS) $(INCLUDES) -o $@ $(word 1,$^)u $(word 2,$^) $(LDFLAGS) $(LIBS)
	rm $(SRCDIR)/$(S_TEST).cu


$(S_TEST)_%: $(S_TEST).c $(LIBDIR)/$(SRC)_%.a
	$(CC) $(CFLAGS) $(INCLUDES) -o $@ $(word 1,$^) $(word 2,$^) $(LDFLAGS) $(LIBS)
custom_test: $(addprefix $(S_TEST)_, $(TYPES)) $(addsuffix $(S_TEST), $(EXTRAS))
	./test_rand_sequential $(n) $(k) 0
	./V1test_rand $(n) $(k) 1
	./V2test_rand $(n) $(k) 1 
	./V3test_rand $(n) $(k) 1
	$(RM) *~ $(addprefix $(S_TEST)_, $(TYPES)) $(addsuffix $(S_TEST), $(EXTRAS)) conf-*

test: $(addprefix $(TEST)_, $(TYPES)) $(addsuffix $(TEST), $(EXTRAS))
	./test_sequential
	./V1test
	./V2test
	./V3test
	$(RM) *~ $(addprefix $(TEST)_, $(TYPES)) $(addsuffix $(TEST), $(EXTRAS))

$(LIBDIR)/$(SRC)_%.a: $(LIBDIR)/$(SRC)_%.o
	ar rcs $@ $<


$(LIBDIR)/$(SRC)_%.o: $(SRCDIR)/$(SRC)_%.c*
	$(NVCC) $(CFLAGS) $(INCLUDES) -o $@ -c $<


clean:
	$(RM) $(SRCDIR)/*.o *~ $(addprefix $(MAIN)_, $(TYPES)) $(addsuffix $(MAIN),$(EXTRAS)) $(addprefix $(TEST)_, $(TYPES)) $(addsuffix $(TEST),$(EXTRAS)) $(LIBDIR)/$(SRC)*.a
