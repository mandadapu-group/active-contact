# adapted from makefile by Prof. Phillip Colella (UC Berkeley)
# 27 Sept. 2018



# compiler flags
CXX = g++
CFLAGS =  -O3
CFLAGS += -std=c++11

# directories for object and dependency files
odir = ./dirO
ddir = ./dirD

# find all source files, generate names of objects and dependencies
SRCFILES := $(notdir $(wildcard ./*.cpp))
OBJS     := $(patsubst %.cpp, $(odir)/%.o, $(SRCFILES))
DEPS     := $(patsubst $(odir)/%.o, $(ddir)/%.d, $(OBJS))

# recipe for object files
$(odir)/%.o:%.cpp GNUmakefile
	mkdir -p $(odir); $(CXX) -c $(CFLAGS) $< -o $@
	mkdir -p $(ddir); $(CXX) -MM $(CFLAGS) $< | sed '1s/^/o.$(DIM)d\//' > $*.d;mv $*.d $(ddir)

-include $(DEPS)

# clean
clean:
	rm -r *.exe $(odir) $(ddir)

# useful tools
listsrc:
	@echo $(SRCFILES)
listobj:
	@echo $(OBJS)
listdep:
	@echo $(DEPS)


# recipe for main method
main: GNUmakefile $(OBJS) 
	$(CXX) $(CFLAGS) $(CALCFLAGS) $(OBJS) -o activeContact.exe


