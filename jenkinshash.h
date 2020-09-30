#ifndef jenkinshash
#define jenkinshash

static inline uint32 JeninksHash(uint32 i)
	{
	i = (i+0x7ed55d16) + (i<<12);
	i = (i^0xc761c23c) ^ (i>>19);
	i = (i+0x165667b1) + (i<<5);
	i = (i+0xd3a2646c) ^ (i<<9);
	i = (i+0xfd7046c5) + (i<<3);
	i = (i^0xb55a4f09) ^ (i>>16);
	return i;
	}

#endif // jenkinshash
