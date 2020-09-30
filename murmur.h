#ifndef murmur_h
#define murmur_h

static inline uint64 murmur64(uint64 h)
	{
	h ^= (h >> 33);
	h *= 0xff51afd7ed558ccdL;
	h ^= (h >> 33);
	h *= 0xc4ceb9fe1a85ec53L;
	h ^= (h >> 33);
	return h;
	}

#endif // murmur_h
