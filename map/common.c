#include "common.h"
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <stddef.h>
#include <dirent.h>
#include <ctype.h>
#include <string.h>
//#include <math.h>



static int numeric_filter (const struct dirent * ent)
{
	int i = 0;
	char c;

	while ((c = ent->d_name[i++]) != '\0') {
		if (! isdigit(c))
			return 0;
	}
	return 1;
}

int get_date_dirs(const char * dir, char ** datedirs[])
{
	struct dirent **eps;
	int numDirs;
	int i;
	char ** dirs;

	numDirs = scandir (dir, &eps, numeric_filter, alphasort);
	if (numDirs > 0) 
	{
		dirs = (char **) malloc(sizeof(char*) * numDirs);
		for (i=0; i<numDirs; ++i) {
			dirs[i] = (char *) malloc(sizeof(char) * strlen(eps[i]->d_name) + 1);
			strcpy(dirs[i], eps[i]->d_name);
		}
		*datedirs = dirs;
	}
	return numDirs;
}

static int flux_filter (const struct dirent * ent)
{
	char * str;

	str = strstr(ent->d_name, "fluxtime");
	if (str != ent->d_name) {
		return 0;
	}
	str = strstr(ent->d_name, ".dat");
	if (str == NULL) {
		return 0;
	}
	if (strlen(str) != 4) {
		return 0;
	}
	return 1;
}

int get_flux_files(const char * dir, char ** fluxfiles[])
{
	struct dirent **eps;
	int num;
	int i;
	char ** files;

	num = scandir (dir, &eps, flux_filter, alphasort);
	if (num > 0) 
	{
		files = (char **) malloc(sizeof(char*) * num);
		for (i=0; i<num; ++i) {
			files[i] = (char *) malloc(sizeof(char) * strlen(eps[i]->d_name) + 1);
			strcpy(files[i], eps[i]->d_name);
		}
		*fluxfiles = files;
	}
	return num;
}
