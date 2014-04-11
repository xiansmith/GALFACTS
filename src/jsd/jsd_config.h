/*
 * Configuration object class.
 *
 * Nov 14, 2006
 * THIS IS FAR FROM COMPLETE
 */


#include <stdio.h>

namespace jsd {

class Config {

	public:
	Config();
	~Config();
	read(char * filename);
	write(char * filename);
	int get_as_int(char * name);
	float get_as_float(char * name);
	char * get_as_string(char * name);

	private:
	
	struct ltstr
	{
		bool operator()(const char* s1, const char* s2) const
		{
			return strcmp(s1, s2) < 0;
		}
	};

	map<const char*, char*, ltstr> pairs;
}

Config::Config()
{
	pairs = new map<const char*, char*, ltstr>();
}

Config::~Config()
{
	delete pairs;
}


char * Config::get_as_string(char * name)
{
	return pairs[name];
}

int Config::get_as_int(char * name)
{
	return atoi(pairs[name]);
}

int Config::get_as_float(char * name)
{
	return strtof(pairs[name], NULL);
}


} //end namespace
