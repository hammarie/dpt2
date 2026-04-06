# src/mpx_dpcr/design.py
import primer3.bindings
import logging

logger = logging.getLogger(__name__)


def design_primers_for_gene(
    seq: str,
    product_size_range=(60, 400),
    num_return=5,
    target=None,
):
    """Design primer pairs plus Primer3-scored metrics."""
    logger.info(f"Designing primers for sequence length: {len(seq)} bp")
    logger.info(f"Product size range: {product_size_range}")

    seq_args = {
        "SEQUENCE_ID": "target",
        "SEQUENCE_TEMPLATE": seq,
    }
    if target:
        seq_args["SEQUENCE_TARGET"] = [target]

    global_args = {
        "PRIMER_TASK": "generic",
        "PRIMER_OPT_SIZE": 20,
        "PRIMER_MIN_SIZE": 18,
        "PRIMER_MAX_SIZE": 25,
        "PRIMER_OPT_TM": 60.0,
        "PRIMER_MIN_TM": 57.0,
        "PRIMER_MAX_TM": 63.0,
        "PRIMER_PAIR_MAX_DIFF_TM": 3.0,
        "PRIMER_MIN_GC": 20.0,
        "PRIMER_MAX_GC": 80.0,
        "PRIMER_MAX_POLY_X": 100,
        "PRIMER_SALT_MONOVALENT": 50.0,
        "PRIMER_DNA_CONC": 50.0,
        "PRIMER_MAX_NS_ACCEPTED": 0,
        "PRIMER_MAX_SELF_ANY": 12,
        "PRIMER_MAX_SELF_END": 8,
        "PRIMER_PAIR_MAX_COMPL_ANY": 12,
        "PRIMER_PAIR_MAX_COMPL_END": 8,
        "PRIMER_PRODUCT_SIZE_RANGE": [product_size_range],
        "PRIMER_NUM_RETURN": num_return,
    }

    logger.debug(f"Primer3 global args: {global_args}")

    results = primer3.bindings.design_primers(seq_args, global_args)

    logger.info(f"Primer3 returned {len(results)} keys")

    # Check for errors
    if "PRIMER_LEFT_0_SEQUENCE" not in results:
        logger.warning(f"No primers designed. Primer3 return keys: {list(results.keys())}")
        # Log all keys to see what we got
        for key in sorted(results.keys()):
            logger.debug(f"  {key}: {results[key]}")
        return []

    # Collect all returned primer pairs
    primer_pairs = []
    for i in range(num_return):
        left_key = f"PRIMER_LEFT_{i}_SEQUENCE"
        right_key = f"PRIMER_RIGHT_{i}_SEQUENCE"
        size_key = f"PRIMER_PAIR_{i}_PRODUCT_SIZE"
        if left_key not in results:
            logger.debug(f"No pair {i}, stopping search")
            break  # stop if fewer pairs were found
        pair = {
            "forward": results[left_key],
            "reverse": results[right_key],
            "product_size": results[size_key],
            "forward_tm": results.get(f"PRIMER_LEFT_{i}_TM"),
            "reverse_tm": results.get(f"PRIMER_RIGHT_{i}_TM"),
            "forward_self_any": results.get(f"PRIMER_LEFT_{i}_SELF_ANY_TH"),
            "forward_self_end": results.get(f"PRIMER_LEFT_{i}_SELF_END_TH"),
            "reverse_self_any": results.get(f"PRIMER_RIGHT_{i}_SELF_ANY_TH"),
            "reverse_self_end": results.get(f"PRIMER_RIGHT_{i}_SELF_END_TH"),
            "pair_compl_any": results.get(f"PRIMER_PAIR_{i}_COMPL_ANY_TH"),
            "pair_compl_end": results.get(f"PRIMER_PAIR_{i}_COMPL_END_TH"),
            "pair_penalty": results.get(f"PRIMER_PAIR_{i}_PENALTY"),
        }
        logger.info(f"Pair {i}: F={pair['forward']} R={pair['reverse']} size={pair['product_size']}")
        primer_pairs.append(pair)

    logger.info(f"Successfully designed {len(primer_pairs)} primer pairs")
    return primer_pairs
