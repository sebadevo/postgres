/*-------------------------------------------------------------------------
 *
 * rangetypes_typanalyze.c
 *	  Functions for gathering statistics from range columns
 *
 * For a range type column, histograms of lower and upper bounds, and
 * the fraction of NULL and empty ranges are collected.
 *
 * Both histograms have the same length, and they are combined into a
 * single array of ranges. This has the same shape as the histogram that
 * std_typanalyze would collect, but the values are different. Each range
 * in the array is a valid range, even though the lower and upper bounds
 * come from different tuples. In theory, the standard scalar selectivity
 * functions could be used with the combined histogram.
 *
 * Portions Copyright (c) 1996-2020, PostgreSQL Global Development Group
 * Portions Copyright (c) 1994, Regents of the University of California
 *
 *
 * IDENTIFICATION
 *	  src/backend/utils/adt/rangetypes_typanalyze.c
 *
 *-------------------------------------------------------------------------
 */
#include "postgres.h"

#include "catalog/pg_operator.h"
#include "commands/vacuum.h"
#include "utils/float.h"
#include "utils/fmgrprotos.h"
#include "utils/lsyscache.h"
#include "utils/rangetypes.h"

static int	float8_qsort_cmp(const void *a1, const void *a2);
static int	range_bound_qsort_cmp(const void *a1, const void *a2, void *arg);
static void compute_range_stats(VacAttrStats *stats,
								AnalyzeAttrFetchFunc fetchfunc, int samplerows, double totalrows);

/*
 * range_typanalyze -- typanalyze function for range columns
 */
Datum
range_typanalyze(PG_FUNCTION_ARGS)
{
	VacAttrStats *stats = (VacAttrStats *) PG_GETARG_POINTER(0);
	TypeCacheEntry *typcache;
	Form_pg_attribute attr = stats->attr;

	/* Get information about range type; note column might be a domain */
	typcache = range_get_typcache(fcinfo, getBaseType(stats->attrtypid));

	if (attr->attstattarget < 0)
		attr->attstattarget = default_statistics_target;

	stats->compute_stats = compute_range_stats;
	stats->extra_data = typcache;
	/* same as in std_typanalyze */
	stats->minrows = 300 * attr->attstattarget;



	PG_RETURN_BOOL(true);
}

/*
 * Comparison function for sorting float8s, used for range lengths.
 */
static int
float8_qsort_cmp(const void *a1, const void *a2)
{
	const float8 *f1 = (const float8 *) a1;
	const float8 *f2 = (const float8 *) a2;

	if (*f1 < *f2)
		return -1;
	else if (*f1 == *f2)
		return 0;
	else
		return 1;
}

/*
 * Comparison function for sorting RangeBounds.
 */
static int
range_bound_qsort_cmp(const void *a1, const void *a2, void *arg)
{
	RangeBound *b1 = (RangeBound *) a1;
	RangeBound *b2 = (RangeBound *) a2;
	TypeCacheEntry *typcache = (TypeCacheEntry *) arg;

	return range_cmp_bounds(typcache, b1, b2);
}

/*
 * compute_range_stats() -- compute statistics for a range column
 */
static void
compute_range_stats(VacAttrStats *stats, AnalyzeAttrFetchFunc fetchfunc,
					int samplerows, double totalrows)
{
	TypeCacheEntry *typcache = (TypeCacheEntry *) stats->extra_data;
	bool		has_subdiff = OidIsValid(typcache->rng_subdiff_finfo.fn_oid);
	int			null_cnt = 0;
	int			non_null_cnt = 0;
	int			non_empty_cnt = 0;
	int			empty_cnt = 0;
	int			range_no;
	int			slot_idx;
	int			num_bins = stats->attr->attstattarget; // 100 pour notre cas
	int			num_hist;

	// Maximum value in the range column, comes from DatumGetInt16
	
	float8	   *lengths,
				*frequencies; // [frequency histogram]
	RangeBound *lowers,
			   *uppers;
	double		total_width = 0;



	/* Allocate memory to hold range bounds and lengths of the sample ranges. */
	lowers = (RangeBound *) palloc(sizeof(RangeBound) * samplerows);
	uppers = (RangeBound *) palloc(sizeof(RangeBound) * samplerows);
	lengths = (float8 *) palloc(sizeof(float8) * samplerows);
	
	fflush(stdout) ;

	/**
	 * [frequency histogram]
	 * We first need to find the min and max bounds of the column
	 */
	RangeBound 	max_bound,
				min_bound;

	// Initial values to intitialize max_bound and min_bound
	RangeType *initRange ;
	
	Datum initValue ;
	bool initIsnull,
		initEmpty ;

	initValue = fetchfunc(stats, 0, &initIsnull);

	initRange = DatumGetRangeTypeP(initValue);

	range_deserialize(typcache, initRange, &min_bound, &max_bound, &initEmpty);

	int16 initialMax = DatumGetInt16(max_bound.val) ;
	int16 initialMin = DatumGetInt16(min_bound.val) ;
	// End of initialization

	/* Loop over the sample ranges. */
	for (range_no = 0; range_no < samplerows; range_no++)
	{
		Datum		value;
		bool		isnull,
					empty;
		RangeType  *range;
		RangeBound	lower,
					upper;
		float8		length,
					frequency; // [frequency histogram]

		vacuum_delay_point();

		value = fetchfunc(stats, range_no, &isnull);
		if (isnull)
		{
			/* range is null, just count that */
			null_cnt++;
			continue;
		}

		/*
		 * XXX: should we ignore wide values, like std_typanalyze does, to
		 * avoid bloating the statistics table?
		 */
		total_width += VARSIZE_ANY(DatumGetPointer(value));


		/* Get range and deserialize it for further analysis. */
		range = DatumGetRangeTypeP(value);
		range_deserialize(typcache, range, &lower, &upper, &empty);

		if (!empty)
		{
			/* Remember bounds and length for further usage in histograms */
			lowers[non_empty_cnt] = lower;
			uppers[non_empty_cnt] = upper;

			/**
			 * [frequency histogram]
			 * Comparison to keep the max bound and min bound
			 */
			if(range_cmp_bounds(typcache, &upper, &max_bound) > 0){
				max_bound = upper ;
			}
			if(range_cmp_bounds(typcache, &lower, &min_bound) < 0){
				min_bound = lower ;
			}

			if (lower.infinite || upper.infinite)
			{
				/* Length of any kind of an infinite range is infinite */
				length = get_float8_infinity();
			}
			else if (has_subdiff)
			{
				/*
				 * For an ordinary range, use subdiff function between upper
				 * and lower bound values.
				 */
				length = DatumGetFloat8(FunctionCall2Coll(&typcache->rng_subdiff_finfo,
														  typcache->rng_collation,
														  upper.val, lower.val));
			}
			else
			{
				/* Use default value of 1.0 if no subdiff is available. */
				length = 1.0;
			}

			// here we fill the length array
			lengths[non_empty_cnt] = length;
			non_empty_cnt++;
		}
		else
			empty_cnt++;

		non_null_cnt++;
	}

	slot_idx = 0;

	/* We can only compute real stats if we found some non-null values. */
	if (non_null_cnt > 0)
	{
		Datum	   *bound_hist_values;
		Datum	   *length_hist_values;
		Datum	   *frequency_hist_values;
		int			pos,
					posfrac,
					delta,
					deltafrac,
					i;
		// Attention à ne pas mettre "double" !!
		float8  	width; // [frequency histogram]

		MemoryContext old_cxt;
		float4	   *emptyfrac;

		stats->stats_valid = true;
		/* Do the simple null-frac and width stats */
		stats->stanullfrac = (double) null_cnt / (double) samplerows;
		stats->stawidth = total_width / (double) non_null_cnt;

		// printf("%d\n", stats->stawidth);
		// fflush(stdout);

		/* Estimate that non-null values are unique */
		stats->stadistinct = -1.0 * (1.0 - stats->stanullfrac);

		/* Must copy the target values into anl_context */
		old_cxt = MemoryContextSwitchTo(stats->anl_context);

		/*
		 * Generate a bounds histogram slot entry if there are at least two
		 * values.
		 */
		if (non_empty_cnt >= 2)
		{
			/* Sort bound values */
			// Le minimum est dans lowers[0] !!
			qsort_arg(lowers, non_empty_cnt, sizeof(RangeBound),
					  range_bound_qsort_cmp, typcache);
			qsort_arg(uppers, non_empty_cnt, sizeof(RangeBound),
					  range_bound_qsort_cmp, typcache);
			fflush(stdout) ;
			num_hist = non_empty_cnt;
			if (num_hist > num_bins)
				num_hist = num_bins + 1;


			bound_hist_values = (Datum *) palloc(num_hist * sizeof(Datum));

			/*
			 * The object of this loop is to construct ranges from first and
			 * last entries in lowers[] and uppers[] along with evenly-spaced
			 * values in between. So the i'th value is a range of lowers[(i *
			 * (nvals - 1)) / (num_hist - 1)] and uppers[(i * (nvals - 1)) /
			 * (num_hist - 1)]. But computing that subscript directly risks
			 * integer overflow when the stats target is more than a couple
			 * thousand.  Instead we add (nvals - 1) / (num_hist - 1) to pos
			 * at each step, tracking the integral and fractional parts of the
			 * sum separately.
			 */
			delta = (non_empty_cnt - 1) / (num_hist - 1);
			deltafrac = (non_empty_cnt - 1) % (num_hist - 1);



			pos = posfrac = 0;

			for (i = 0; i < num_hist; i++)
			{

				// Bound hist content is a pointer to a range
				bound_hist_values[i] = PointerGetDatum(range_serialize(typcache,
																	   &lowers[pos],
																	   &uppers[pos],
																	   false));

				pos += delta;
				posfrac += deltafrac;
				if (posfrac >= (num_hist - 1))
				{
					/* fractional part exceeds 1, carry to integer part */
					pos++;
					posfrac -= (num_hist - 1);
				}
			}
			


			stats->stakind[slot_idx] = STATISTIC_KIND_BOUNDS_HISTOGRAM;
			stats->stavalues[slot_idx] = bound_hist_values;

			/*
			printf("%d \n", *bound_hist_values) ;
			fflush(stdout) ;
			*/

			stats->numvalues[slot_idx] = num_hist;
			slot_idx++;
		}

		/*
		 * Generate a length histogram slot entry if there are at least two
		 * values.
		 */
		if (non_empty_cnt >= 2)
		{
			/*
			 * Ascending sort of range lengths for further filling of
			 * histogram
			 */
			qsort(lengths, non_empty_cnt, sizeof(float8), float8_qsort_cmp);

			num_hist = non_empty_cnt;
			if (num_hist > num_bins)
				num_hist = num_bins + 1;

			length_hist_values = (Datum *) palloc(num_hist * sizeof(Datum));

			/*
			 * The object of this loop is to copy the first and last lengths[]
			 * entries along with evenly-spaced values in between. So the i'th
			 * value is lengths[(i * (nvals - 1)) / (num_hist - 1)]. But
			 * computing that subscript directly risks integer overflow when
			 * the stats target is more than a couple thousand.  Instead we
			 * add (nvals - 1) / (num_hist - 1) to pos at each step, tracking
			 * the integral and fractional parts of the sum separately.
			 */
			delta = (non_empty_cnt - 1) / (num_hist - 1);
			// pour un tableau de 100 000 lignes, num_hist sera 101 et alors delta entier sera 999
			deltafrac = (non_empty_cnt - 1) % (num_hist - 1);
			// la partie décimale de delta sera de 99
			
			pos = posfrac = 0;

			for (i = 0; i < num_hist; i++)
			{
				length_hist_values[i] = Float8GetDatum(lengths[pos]);

				pos += delta;
				posfrac += deltafrac;
				if (posfrac >= (num_hist - 1))
				{
					/* fractional part exceeds 1, carry to integer part */
					pos++;
					posfrac -= (num_hist - 1);
				}
			}
			fflush(stdout) ;
			
		}
		else
		{
			/*
			 * Even when we don't create the histogram, store an empty array
			 * to mean "no histogram". We can't just leave stavalues NULL,
			 * because get_attstatsslot() errors if you ask for stavalues, and
			 * it's NULL. We'll still store the empty fraction in stanumbers.
			 */
			length_hist_values = palloc(0);
			num_hist = 0;
		}
		stats->staop[slot_idx] = Float8LessOperator;
		stats->stacoll[slot_idx] = InvalidOid;
		stats->stavalues[slot_idx] = length_hist_values;
		stats->numvalues[slot_idx] = num_hist;
		stats->statypid[slot_idx] = FLOAT8OID;
		stats->statyplen[slot_idx] = sizeof(float8);
		stats->statypbyval[slot_idx] = FLOAT8PASSBYVAL;
		stats->statypalign[slot_idx] = 'd';

		/* Store the fraction of empty ranges */
		emptyfrac = (float4 *) palloc(sizeof(float4));
		*emptyfrac = ((double) empty_cnt) / ((double) non_null_cnt);
		stats->stanumbers[slot_idx] = emptyfrac;
		stats->numnumbers[slot_idx] = 1;

		stats->stakind[slot_idx] = STATISTIC_KIND_RANGE_LENGTH_HISTOGRAM;
		slot_idx++;


		/*
		 * Generate a [frequency histogram] slot entry if there are at least two
		 * values.
		 * 
		 * The frequency histogram carries two informations :
		 * - The width of each bin : stored at index 0
		 * - The values of the frequency for each bin : stored in the indexes 1 to (num_hist+1)
		 */
		if (non_empty_cnt >= 2)
		{

			// Similarly to length histogram and bounds histograms
			num_hist = non_empty_cnt;
			if (num_hist > num_bins)
				num_hist = num_bins + 1;
			
			
			/**
			 * num_hist bins to store the frequencies, 1 bin to store the width
			 * Frequencies are stored as floats
			 */
			frequency_hist_values = (Datum *) palloc((num_hist+1) * sizeof(Datum));
			frequencies = palloc(num_hist*sizeof(float8)) ;


			float8 frequency_at_index,
					column_length ;
			
			// FH initialization
			for (int i = 0; i < num_hist; i++)
			{
				frequencies[i] = 0 ;
			}
			
			// We compute length by comparing max_bound and min_bound, previously stored
			column_length = DatumGetFloat8(FunctionCall2Coll(&typcache->rng_subdiff_finfo,
														  typcache->rng_collation,
														  max_bound.val, min_bound.val)); 
			
			// we take the entire part of the division
			width = div(column_length,num_hist).quot  ;
			// this avoids going out of bounds, even with reminder is 0
			width ++ ;

			/**
			 * Need to iterate over the ranges of the sample rows
			 * to locate the bins of the FH to increment
			 */
			for (int range_no = 0; range_no < samplerows ; range_no++)
			{
				
				RangeType *current_range ;
				RangeBound	current_lower,
							current_upper;
				Datum current_value ;
				bool current_isnull,
					current_isEmpty ;

				current_value = fetchfunc(stats, range_no, &initIsnull);
				current_range = DatumGetRangeTypeP(current_value);
				range_deserialize(typcache, current_range, &current_lower, &current_upper, &current_isEmpty);


				// Operations between the current range and the max_bound to find the relative position
				float8 compLower, compUpper ;
				compLower = DatumGetFloat8(FunctionCall2Coll(&typcache->rng_subdiff_finfo,
														  typcache->rng_collation,
														  max_bound.val, current_lower.val)); 
				compUpper = DatumGetFloat8(FunctionCall2Coll(&typcache->rng_subdiff_finfo,
														  typcache->rng_collation,
														  max_bound.val, current_upper.val)); 


				// Integer division to set bound to increment the histogram
				int delta_low = div((column_length-compLower), width).quot ;
				div_t delta_up_div = div((column_length-compUpper), width) ;
				int delta_up = delta_up_div.quot ;
				int rem = delta_up_div.rem; 

				if(rem != 0){
					delta_up ++ ;
				}

				// Increments the FH for each bin in the covered bins
				for (int i = delta_low ; i < delta_up ; i++)
				{
					// Frequencies array is an array of float, we can just increment the values
					frequencies[i]++ ;
				}
				
			}

			// Store the width pointer in the first element of the histogram
			frequency_hist_values[0] = Float8GetDatum(width) ;

			/**
			 * Index = 0 has already been treated. For index > 0 we store the pointers
			 * to the float frequencies in the frequencies array
			 */
			for (int idx = 1; idx < num_hist+1; idx++)
			{
				frequency_hist_values[idx] = Float8GetDatum(frequencies[idx-1]) ;
			}

			
			 // Print of frequency histogram
			
			printf("frequence_histogram = WIDTH %f \n[", width);
			for (i = 0; i < num_hist+1; i++)
			{
				int frequency = DatumGetFloat8(frequency_hist_values[i]) ;
				printf("%d", frequency) ;
				if (i < num_hist+1 - 1)
					printf(", ");
			}
			printf("]\n") ;
			
			fflush(stdout) ;
		}
		else
		{
			/*
			 * Even when we don't create the histogram, store an empty array
			 * to mean "no histogram". We can't just leave stavalues NULL,
			 * because get_attstatsslot() errors if you ask for stavalues, and
			 * it's NULL. We'll still store the empty fraction in stanumbers.
			 */
			frequency_hist_values = palloc(0);
			num_hist = 0;
		}

		

		stats->staop[slot_idx] = Float8LessOperator;
		stats->stacoll[slot_idx] = InvalidOid;
		stats->stavalues[slot_idx] = frequency_hist_values;
		stats->numvalues[slot_idx] = num_hist+1; // we add 1 for the width value at the index 0
		stats->statypid[slot_idx] = FLOAT8OID;
		stats->statyplen[slot_idx] = sizeof(float8);
		stats->statypbyval[slot_idx] = FLOAT8PASSBYVAL;
		stats->statypalign[slot_idx] = 'd';

		/* Store the fraction of empty ranges */
		emptyfrac = (float4 *) palloc(sizeof(float4));
		*emptyfrac = ((double) empty_cnt) / ((double) non_null_cnt);
		stats->stanumbers[slot_idx] = emptyfrac;
		stats->numnumbers[slot_idx] = 1;

		stats->stakind[slot_idx] = STATISTIC_KIND_FREQUENCY_HISTOGRAM;
		slot_idx++;

		

		MemoryContextSwitchTo(old_cxt);


		
	}
	else if (null_cnt > 0)
	{
		/* We found only nulls; assume the column is entirely null */
		stats->stats_valid = true;
		stats->stanullfrac = 1.0;
		stats->stawidth = 0;	/* "unknown" */
		stats->stadistinct = 0.0;	/* "unknown" */
	}

	/*
	 * We don't need to bother cleaning up any of our temporary palloc's. The
	 * hashtable should also go away, as it used a child memory context.
	 */



}