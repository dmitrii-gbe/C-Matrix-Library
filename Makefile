ifeq ($(origin CC),default)
	CC = gcc
endif

CFLAGS ?= -O3 -Wall -Werror -Wextra -std=c11
CFLAGSDEBUG ?= -g -Wall -Werror -Wextra -std=c11
# CFLAGS ?= -g -std=c11
COBJ_DIR = ./object_files
BINDIR := .
STATIC_LIB = s21_matrix.a
TEST = ./test
LIBS = -lcheck
CSRC = s21_matrix.c
TEST_CSRC := $(TEST)/test.c
COBJ := $(addprefix $(COBJ_DIR)/,$(CSRC:.c=.o))
TEST_COBJ := $(TEST_CSRC:.c=.o)
DEPS := $(COBJ:.o=.d)
TEST_DEPS := $(TEST_COBJ:.o=.d)

ifeq ($(shell uname -s),Linux)
	LIBS += -lsubunit
endif

ifeq ($(shell uname -s),Linux)
	LCOV_DIR = ${TEST}
	OPENER := xdg-open
else ifeq ($(shell uname -s),Darwin)
	LCOV_DIR = ${BINDIR}
	OPENER := open
endif

.PHONY: all clean test gcov_report clang-format

all: s21_matrix.a test gcov_report

$(COBJ) : $(COBJ_DIR)/%.o : %.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(TEST_COBJ) : $(TEST)/%.o : $(TEST)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(DEPS) : $(COBJ_DIR)/%.d : %.c
	@mkdir -p $(@D)
	$(CC) -E $(CFLAGS) $< -MM -MT $(@:.d=.o) > $@

$(TEST_DEPS) : $(TEST)/%.d : $(TEST)/%.c
	@mkdir -p $(@D)
	$(CC) -E $(CFLAGS) $< -MM -MT $(@:.d=.o) > $@

$(BINDIR)/$(STATIC_LIB): $(COBJ) 
	ar -rcs $@ $^ 

test: $(TEST_COBJ) $(BINDIR)/$(STATIC_LIB)
	$(CC) $(CFLAGS) $< -L. $(STATIC_LIB) $(LIBS) -o $(TEST)/test
	@chmod +x ${TEST}/test
	@cd ${TEST} && ./test 

gcov_report: $(BINDIR)/$(CSRC) $(TEST_CSRC)
	@make clean
	@$(CC) -fprofile-arcs -ftest-coverage $(CFLAGSDEBUG) $^ -o $(TEST)/test $(LIBS)
	@chmod +x ${TEST}/test
	@cd ${TEST} && ./test 
	lcov -t "gcov_report" -o gcov_report.info -c -d ${LCOV_DIR}
	genhtml -o gcov_report gcov_report.info 
	${OPENER} gcov_report/index.html


clean:
	@rm -rf $(COBJ_DIR)
	@rm -f ./*.a
	@rm -f test/*.gcda
	@rm -f test/*.gcno
	@rm -f *.gcda
	@rm -f *.gcno
	@rm -f *.info
	@rm -rf test/gcov_report/
	@rm -rf gcov_report/
	@rm -f $(TEST)/test $(TEST)/test.o $(TEST)/test.d

clang-format:
	clang-format -i *.c *.h test/*.c test/*.h

NODEPS = clean clang-format

ifeq (0, $(words $(findstring $(MAKECMDGOALS), $(NODEPS))))
	include $(TEST_DEPS)
	include $(DEPS)
endif

rebuild:
	@make clean
	@make all
