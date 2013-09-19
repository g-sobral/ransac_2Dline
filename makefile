##############################################################################
SHELL = /bin/sh

CC := gcc
CFLAGS := -std=c99 -Wall -pedantic -O3
LDFLAGS := -lm

SRCDIR := src
BUILDDIR := build
SRCEXT := c
EXEC := ransac2D

SOURCES := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.o))

all : $(EXEC)

$(EXEC): $(OBJECTS)
	@echo " Linking..."; $(CC) -o $@ $^ $(LDFLAGS)

$(BUILDDIR)/%.o: $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(BUILDDIR)
	@echo " CC $<"; $(CC) $(CFLAGS) -c -o $@ $<

clean:
	@echo " Cleaning..."; rm -r $(BUILDDIR) $(EXEC)

##############################################################################
