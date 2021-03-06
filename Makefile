### CONFIGURE - EDIT TO PREFERENCE

# Compiler, flags and include directories (-I)
CC    := g++
CCFL  := -Wall -Werror -pedantic -std=c++11 -O3
INC   := lib/boost_1_67_0

# Linker, flags and external libraries (-l) 
CCLD  := $(CC)
CCLDF := 
LIB   := z

# Executable name, source dir and compile dir
EXEC  := qgs
SRC   := src
BIN   := bin

# Colours used in messages
C_DEFAULT := \033[m
C_CORRECT := \033[0;32m
C_FTERROR := \033[0;31m
C_WARNING := \033[0;33m
C_MESSAGE := \033[0;34m

# Toggle fancy printing (colours and formatting) [true/false]
FANCY = false

### CANNED RECIPIES FOR COMPILING AND LINKING - EDIT TO PREFERENCE

# Recipe for compiling
ifeq ($(FANCY),true)
define compile
	@echo -en "$(C_MESSAGE)----- [$(N)/$(words $(OBJ))] $< ----- $(C_DEFAULT)"
	@$(CC) -c $< -o $@ $(CCFL) $(INC) 2> .temp.log || touch .temp.err
	$(result_eval)
endef
else
define compile
	@echo -n "[$(N)/$(words $(OBJ))] "
	$(CC) -c $< -o $@ $(CCFL) $(INC)
endef
endif

# Recipe for linking
ifeq ($(FANCY),true)
define link
	@echo -ne "$(C_MESSAGE)----- LINKING ----- $(C_DEFAULT)"
	@$(CCLD) $(OBJ) -o $(EXEC) $(CCLDF) $(LIB) 2> .temp.log || touch .temp.err
	$(result_eval)
endef
else
define link
	@echo -n "[LINKING] "
	$(CCLD) $(OBJ) -o $(EXEC) $(CCLDF) $(LIB)
endef
endif

# Recipe for result evaluation, used for fancy printing
define result_eval
	@if test -e .temp.err; then echo -e "$(C_FTERROR)[ERR]"; elif test -s .temp.log; then echo -e "$(C_WARNING)[WARN]"; else echo -e "$(C_CORRECT)[OK]"; fi;
	@cat .temp.log && echo -en "$(C_DEFAULT)"
	@if test -e .temp.err; then rm -f .temp.log .temp.err && false; fi;
	@rm -f .temp.log
endef

### DO NOT CHANGE ANYTHING FROM THIS POINT ON, INCLUDING THIS LINE AND OTHER COMMENTS

### RULES FOR TRANSFORMING COMPILER ARGUMENTS

# Adding -I to $(INC)
INC := $(patsubst %,-I%,$(INC))

# Adding -l to $(LIB)
LIB := $(patsubst %,-l%,$(LIB))

### OBJECT FILE DATA

# Object files
OBJ   := $(BIN)/boost_iostreams/gzip.o $(BIN)/boost_iostreams/zlib.o $(BIN)/command_line_arguments.o $(BIN)/gene_score.o $(BIN)/plinkbedreader.o $(BIN)/plinkdosagereader.o $(BIN)/qgs.o $(BIN)/snpreader.o $(BIN)/vcfreader.o

### PHONY RULES

.PHONY	: all clean createbin

# Main goal
all : createbin $(OBJ)
	$(link)

# Clean executable and object file directories
clean :
	rm -fr $(EXEC) $(BIN)

# Create object file directories
createbin :
	@mkdir -p $(sort $(dir $(OBJ)))

### OBJECT FILE RECIPIES FOR COMPILING

$(BIN)/boost_iostreams/gzip.o : $(SRC)/boost_iostreams/gzip.cpp 
	$(eval N := 1)
	$(compile)

$(BIN)/boost_iostreams/zlib.o : $(SRC)/boost_iostreams/zlib.cpp 
	$(eval N := 2)
	$(compile)

$(BIN)/command_line_arguments.o : $(SRC)/command_line_arguments.cc src/cmd_line_options.h src/command_line_arguments.h src/gplv3.h src/gzfile.h src/log.h
	$(eval N := 3)
	$(compile)

$(BIN)/gene_score.o : $(SRC)/gene_score.cc src/gene_score.h src/gzfile.h src/log.h src/snpreader.h
	$(eval N := 4)
	$(compile)

$(BIN)/plinkbedreader.o : $(SRC)/plinkbedreader.cc src/gzfile.h src/log.h src/plinkbedreader.h src/snpreader.h
	$(eval N := 5)
	$(compile)

$(BIN)/plinkdosagereader.o : $(SRC)/plinkdosagereader.cc src/gzfile.h src/log.h src/naturalsort.h src/plinkdosagereader.h src/snpreader.h
	$(eval N := 6)
	$(compile)

$(BIN)/qgs.o : $(SRC)/qgs.cc src/command_line_arguments.h src/genblock.h src/gene_score.h src/gzfile.h src/log.h src/plinkbedreader.h src/plinkdosagereader.h src/qgs.h src/snpreader.h src/vcfreader.h
	$(eval N := 7)
	$(compile)

$(BIN)/snpreader.o : $(SRC)/snpreader.cc src/gzfile.h src/log.h src/snpreader.h
	$(eval N := 8)
	$(compile)

$(BIN)/vcfreader.o : $(SRC)/vcfreader.cc src/gzfile.h src/log.h src/snpreader.h src/vcfreader.h
	$(eval N := 9)
	$(compile)

