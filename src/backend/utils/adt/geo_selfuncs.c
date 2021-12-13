/*-------------------------------------------------------------------------
 *
 * geo_selfuncs.c
 *	  Selectivity routines registered in the operator catalog in the
 *	  "oprrest" and "oprjoin" attributes.
 *
 * Portions Copyright (c) 1996-2020, PostgreSQL Global Development Group
 * Portions Copyright (c) 1994, Regents of the University of California
 *
 *
 * IDENTIFICATION
 *	  src/backend/utils/adt/geo_selfuncs.c
 *
 *	XXX These are totally bogus.  Perhaps someone will make them do
 *	something reasonable, someday.
 *
 *-------------------------------------------------------------------------
 */
#include "postgres.h"

#include "utils/builtins.h"
#include "utils/geo_decls.h"
#include "access/htup_details.h"
#include "catalog/pg_statistic.h"
#include "nodes/pg_list.h"
#include "optimizer/pathnode.h"
#include "optimizer/optimizer.h"
#include "utils/lsyscache.h"
#include "utils/typcache.h"
#include "utils/selfuncs.h"
#include "utils/rangetypes.h"

/*
 *	Selectivity functions for geometric operators.  These are bogus -- unless
 *	we know the actual key distribution in the index, we can't make a good
 *	prediction of the selectivity of these operators.
 *
 *	Note: the values used here may look unreasonably small.  Perhaps they
 *	are.  For now, we want to make sure that the optimizer will make use
 *	of a geometric index if one is available, so the selectivity had better
 *	be fairly small.
 *
 *	In general, GiST needs to search multiple subtrees in order to guarantee
 *	that all occurrences of the same key have been found.  Because of this,
 *	the estimated cost for scanning the index ought to be higher than the
 *	output selectivity would indicate.  gistcostestimate(), over in selfuncs.c,
 *	ought to be adjusted accordingly --- but until we can generate somewhat
 *	realistic numbers here, it hardly matters...
 */


/*
 * Selectivity for operators that depend on area, such as "overlap" (well not anymore).
 */

Datum
areasel(PG_FUNCTION_ARGS)
{
	PG_RETURN_FLOAT8(0.005);
}

Datum
areajoinsel(PG_FUNCTION_ARGS)
{
	PG_RETURN_FLOAT8(0.005);
}

/*
 *	positionsel
 *
 * How likely is a box to be strictly left of (right of, above, below)
 * a given box?
 */

Datum
positionsel(PG_FUNCTION_ARGS)
{
	PG_RETURN_FLOAT8(0.1);
}

Datum
positionjoinsel(PG_FUNCTION_ARGS)
{
	PG_RETURN_FLOAT8(0.1);
}

/*
 *	contsel -- How likely is a box to contain (be contained by) a given box?
 *
 * This is a tighter constraint than "overlap", so produce a smaller
 * estimate than areasel does.
 */

Datum
contsel(PG_FUNCTION_ARGS)
{
	PG_RETURN_FLOAT8(0.001);
}

Datum
contjoinsel(PG_FUNCTION_ARGS)
{
	PG_RETURN_FLOAT8(0.001);
}

float8 join_estimation(int nhist_small, int nhist_big, AttStatsSlot sslot_freq_small, AttStatsSlot sslot_freq_big, float8 small_width, float8 biggest_width, float8 distance) ;

/*
 * Range Overlaps Join Selectivity.
 */
Datum
rangeoverlapsjoinsel(PG_FUNCTION_ARGS)
{
    PlannerInfo *root = (PlannerInfo *) PG_GETARG_POINTER(0);
    Oid         operator = PG_GETARG_OID(1);
    List       *args = (List *) PG_GETARG_POINTER(2);
    JoinType    jointype = (JoinType) PG_GETARG_INT16(3);
    SpecialJoinInfo *sjinfo = (SpecialJoinInfo *) PG_GETARG_POINTER(4);
    Oid         collation = PG_GET_COLLATION();

    double      selec = 0.005;
    float8      cardinality_estimation = 0; // number of rows that will be estimated.

    VariableStatData vardata1, 
                vardata2; 


    Oid         opfuncoid;

    // We need : the frequency histograms of the two columns and their min range bound.
    AttStatsSlot sslot_freq1, 
                    sslot_freq2,
                    sslot_bound1,
                    sslot_bound2;

    
    // Size of the frequency histograms : also contains the width => to decrement of 1 when using.
    int         nhist1, nhist2;
    // Bounds that will be compared to know which one is the smallest (important).
    RangeBound min_bound1, min_bound2 ;
    // Useless bounds that need to be here for deserialization.
    RangeBound upper_bound1, upper_bound2 ;
    


    // column with the smallest min range bound. 
    Datum     *frequency_hist_smallest_min;
    Datum     *frequency_hist_biggest_min;
    int         i;
    Form_pg_statistic stats_freq1 = NULL,
                        stats_freq2 = NULL ,
                        stats_bound1 = NULL ,
                        stats_bound2 = NULL ;

    TypeCacheEntry *typcache = NULL;
    bool        join_is_reversed;
    bool        empty;


    get_join_variables(root, args, sjinfo,
                       &vardata1, &vardata2, &join_is_reversed);

    typcache = range_get_typcache(fcinfo, vardata1.vartype);
    opfuncoid = get_opcode(operator);

    memset(&sslot_freq1, 0, sizeof(sslot_freq1));
    memset(&sslot_freq2, 0, sizeof(sslot_freq2));
    memset(&sslot_bound1, 0, sizeof(sslot_bound1));
    memset(&sslot_bound2, 0, sizeof(sslot_bound2));
    

    /* Can't use the histogram with insecure range support functions. */
    if (!statistic_proc_security_check(&vardata1, opfuncoid))
        PG_RETURN_FLOAT8((float8) selec);

    if (HeapTupleIsValid(vardata1.statsTuple))
    {
        stats_freq1 = (Form_pg_statistic) GETSTRUCT(vardata1.statsTuple);
        stats_bound1 = (Form_pg_statistic) GETSTRUCT(vardata1.statsTuple);

        /* Try to get fraction of empty ranges. */
        if (!get_attstatsslot(&sslot_freq1, vardata1.statsTuple,
                             STATISTIC_KIND_FREQUENCY_HISTOGRAM,
                             InvalidOid, ATTSTATSSLOT_VALUES))
        {
            ReleaseVariableStats(vardata1);
            PG_RETURN_FLOAT8((float8) selec);
        }
        if (!get_attstatsslot(&sslot_bound1, vardata1.statsTuple,
                             STATISTIC_KIND_BOUNDS_HISTOGRAM,
                             InvalidOid, ATTSTATSSLOT_VALUES))
        {
            ReleaseVariableStats(vardata1);
            PG_RETURN_FLOAT8((float8) selec);
        }
        
    }

    if(HeapTupleIsValid(vardata2.statsTuple))
    {
        stats_freq2 = (Form_pg_statistic) GETSTRUCT(vardata2.statsTuple);
        stats_bound2 = (Form_pg_statistic) GETSTRUCT(vardata2.statsTuple);

        if (!get_attstatsslot(&sslot_freq2, vardata2.statsTuple,
                             STATISTIC_KIND_FREQUENCY_HISTOGRAM,
                             InvalidOid, ATTSTATSSLOT_VALUES))
        {
            ReleaseVariableStats(vardata2);
            PG_RETURN_FLOAT8((float8) selec);
        }
        if (!get_attstatsslot(&sslot_bound2, vardata2.statsTuple,
                             STATISTIC_KIND_BOUNDS_HISTOGRAM,
                             InvalidOid, ATTSTATSSLOT_VALUES))
        {
            ReleaseVariableStats(vardata2);
            PG_RETURN_FLOAT8((float8) selec);
        }
    }

    // Get size of the histograms (careful to decrement 1 when needed !).
    nhist1 = sslot_freq1.nvalues;
    nhist2 = sslot_freq2.nvalues;
    
    // Cartesian product cardinality
    double rows1 = vardata1.rel->rows ;
    double rows2 = vardata2.rel->rows ;
    double cartesian_product_cardinality = rows1*rows2 ;

    // Gets the min_bound of the first table.
    range_deserialize(typcache, DatumGetRangeTypeP(sslot_bound1.values[0]),
        &min_bound1, &upper_bound1, &empty) ;
    if (empty)
        elog(ERROR, "bounds histogram contains an empty range");

    // Gets the min_bound of the second table.
        range_deserialize(typcache, DatumGetRangeTypeP(sslot_bound2.values[0]),
        &min_bound2, &upper_bound2, &empty) ;
    if (empty)
        elog(ERROR, "bounds histogram contains an empty range");

    // Gets the rows width of the frequency histograms.
    float8 width_1 = DatumGetFloat8(sslot_freq1.values[0]);
    float8 width_2 = DatumGetFloat8(sslot_freq2.values[0]);
    
    // computes the distance between the minimum of each tables.
    float8 distance = abs(DatumGetFloat8(FunctionCall2Coll(&typcache->rng_subdiff_finfo,
														  typcache->rng_collation,
														  min_bound1.val, min_bound2.val)));

    // We try to find which tables has the lowestt minimum
    int cmp_res = range_cmp_bounds(typcache, &min_bound1, &min_bound2) ;
    if( cmp_res < 0){
        // minbound 1 < minbound2, so table 1 has the minimum value
        cardinality_estimation = join_estimation(nhist1, nhist2, sslot_freq1, sslot_freq2, width_1, width_2, distance) ;
    }
    else{
        // minbound 2 <= minbound1
        cardinality_estimation = join_estimation(nhist2, nhist1, sslot_freq2, sslot_freq1, width_2, width_1, distance) ;
    }

    selec = (Selectivity) cardinality_estimation/(Selectivity)cartesian_product_cardinality ;
    free_attstatsslot(&sslot_freq1);
    free_attstatsslot(&sslot_freq2);
    free_attstatsslot(&sslot_bound1);
    free_attstatsslot(&sslot_bound2);

    ReleaseVariableStats(vardata1);
    ReleaseVariableStats(vardata2);

    CLAMP_PROBABILITY(selec);
    PG_RETURN_FLOAT8((float8) selec);
}

/*
 * Computes the join estimation value for range overlaps "&&".
 */
float8 join_estimation(int nhist_small, int nhist_big, AttStatsSlot sslot_freq_small, AttStatsSlot sslot_freq_big, float8 small_width, float8 biggest_width, float8 distance){
    float8 cardinality_estimation = 0 ;
    Datum *frequency_hist_smallest_min, *frequency_hist_biggest_min ;
    frequency_hist_smallest_min = (Datum *) palloc(sizeof(Datum) * (nhist_small-1));
    frequency_hist_biggest_min = (Datum *) palloc(sizeof(Datum) * (nhist_big-1));

    // Histogram filling
    for (int i = 0 ; i < nhist_small-1 ; i ++){
        frequency_hist_smallest_min[i] = sslot_freq_small.values[i+1] ;
    }
    for (int i = 0 ; i < nhist_big-1 ; i ++){
        frequency_hist_biggest_min[i] = sslot_freq_big.values[i+1] ;
    }

    /*
     * To compute the value of the estimation we will go trough the histogram with the smallest value.
     * and for each slot of the frequency histogram we will compute wich slots it crosses with the second frequency histogram.
     * to do so, we determine de first index 'beg' and the last index it croses 'end'.
     * We will then multiply each rows of the second frequency histogram in the interval [beg, end] by the current value of the first
     * frequency histogram. This result will be an estimation (usualy much bigger than the actual join) of the numbers of rows of the join
     * table.
     */
    for (int i = 0; i < nhist_small-1; i++)
    {
        int beg = (i*small_width-distance)/biggest_width ;

        div_t div_end = div((i+1)*small_width-distance, biggest_width) ;
        int end = div_end.quot  ;
        if(!div_end.rem ){
            end ++ ;
        }

        for (int j = beg; j < end; j++)
        {
            if (j > -1 && j < nhist_big-1)
            {
                double a = DatumGetFloat8(frequency_hist_smallest_min[i]) ;
                double b = DatumGetFloat8(frequency_hist_biggest_min[j]) ;
                double c = a*b ;
                cardinality_estimation += c ;
            }
        }
    }


    double total_sum = 0 ;

    // gets the theoretical total number of rows of the frequency histogram to compute the average numbers of rows per bins
    // in the frequency histogram.
    for (int i = 0; i < nhist_small; i++)
    {
        total_sum += DatumGetFloat8(frequency_hist_smallest_min[i]) ;
    }

    //computes the entropy of the histogram with the lowest value and the theoretical maximum entropy of that histogram 
    float8 entropy_low = 0 ;
	float8 max_entropy_low = 0 ;
	float8 probability_low = (float8) 1/nhist_small ;
    for (int i = 0; i < nhist_small; i++)
	{
		float8 current_frequency = DatumGetFloat8(frequency_hist_smallest_min[i]) ;
		float8 current_probability = current_frequency/total_sum;

		max_entropy_low -= probability_low*log2(probability_low) ;

		if(current_frequency>0){
			entropy_low -= current_probability*log2(current_probability) ;
		}
	}


    total_sum = 0 ; 

    for (int i = 0; i < nhist_big; i++)
    {
        total_sum += DatumGetFloat8(frequency_hist_biggest_min[i]) ;
    }

    // Entropy
	float8 entropy_high = 0 ;
	float8 max_entropy_high = 0 ;
	float8 probability_high = (float8) 1/nhist_big ;


	for (int i = 0; i < nhist_big; i++)
	{
		float8 current_frequency = DatumGetFloat8(frequency_hist_biggest_min[i]) ;
		float8 current_probability = current_frequency/ total_sum ;

		max_entropy_high -= probability_high*log2(probability_high) ;

		if(current_frequency>0){
			entropy_high -= current_probability*log2(current_probability) ;
		}
	}

    double dampening_factor = pow((max_entropy_high/entropy_high) * (max_entropy_low/entropy_low), 4.25) ;

    return cardinality_estimation/dampening_factor ;
}