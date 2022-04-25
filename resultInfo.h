#pragma once
#include  <memory>

class ResultInfo
{
private:
	double __RunningTime = -1.0, __Influence = -1.0, __InfluenceOriginal = -1.0, __Approx = -1.0;
	int __SeedSize = 0;
	size_t __RRsetsSize = 0;
	Nodelist __VecSeed;
public:
	ResultInfo()
	{
	}

	~ResultInfo()
	{
	}

	/// Get running time.
	double get_running_time() const
	{
		return __RunningTime;
	}

	/// Get influence spread.
	double get_influence() const
	{
		return __Influence;
	}

	/// Get self-estimated influence spread.
	double get_influence_original() const
	{
		return __InfluenceOriginal;
	}

	/// Get approximation guarantee.
	double get_approximation() const
	{
		return __Approx;
	}

	/// Get seed sets.
	Nodelist get_seed_vec() const
	{
		return __VecSeed;
	}

	/// Get seed size.
	int get_seed_size() const
	{
		return __SeedSize;
	}

	/// Get the number of RR sets.
	size_t get_RRsets_size() const
	{
		return __RRsetsSize;
	}

	/// Set running time.
	void set_running_time(const double value)
	{
		__RunningTime = value;
	}

	/// Set influence spread.
	void set_influence(const double value)
	{
		__Influence = value;
	}

	/// Set self-estimated influence spread
	void set_influence_original(const double value)
	{
		__InfluenceOriginal = value;
	}

	/// Set approximation guarantee.
	void set_approximation(const double value)
	{
		__Approx = value;
	}

	/// Set seed sets.
	void set_seed_vec(Nodelist& vecSeed)
	{
		__VecSeed = vecSeed;
		set_seed_size((int)vecSeed.size());
	}

	/// Set seed size.
	void set_seed_size(const int value)
	{
		__SeedSize = value;
	}

	/// Set the number of RR sets.
	void set_RR_sets_size(const size_t value)
	{
		__RRsetsSize = value;
	}
};

typedef ResultInfo TResult;
typedef std::shared_ptr<ResultInfo> PResult;
