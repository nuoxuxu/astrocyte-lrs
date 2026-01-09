# unctions with no imports within SQANTI3
import bisect
from collections.abc import Iterable
import math
import os
import re
import polars as pl

from contextlib import contextmanager
@contextmanager
def Timer(message=None, verbose=True):
    """
    Use "with Timer(message):" to time the code inside the with block. Based on
    preshing.com/20110924/timing-your-code-using-pythons-with-statement
    
    Args:
        message: a message to print when starting the with block (with "..."
                 after) and ending the with block (with the time after)
        verbose: if False, disables the Timer. This is useful to conditionally
                 run the Timer based on the value of a boolean variable.
    """
    if verbose:
        from timeit import default_timer
        if message is not None:
            print(f'{message}...')
        start = default_timer()
        aborted = False
        try:
            yield
        except Exception as e:
            aborted = True
            raise e
        finally:
            end = default_timer()
            duration = end - start
            
            days = int(duration // 86400)
            hours = int((duration % 86400) // 3600)
            minutes = int((duration % 3600) // 60)
            seconds = int(duration % 60)
            milliseconds = int((duration * 1000) % 1000)
            microseconds = int((duration * 1000000) % 1000)
            nanoseconds = int((duration * 1000000000) % 1000)
            
            time_parts = []
            if days > 0:
                time_parts.append(f'{days} {plural("day", days)}')
            if hours > 0:
                time_parts.append(f'{hours}h')
            if minutes > 0:
                time_parts.append(f'{minutes}m')
            if seconds > 0:
                time_parts.append(f'{seconds}s')
            if milliseconds > 0:
                time_parts.append(f'{milliseconds}ms')
            if microseconds > 0:
                time_parts.append(f'{microseconds}µs')
            if nanoseconds > 0:
                time_parts.append(f'{nanoseconds}ns')
            
            time_str = \
                ' '.join(time_parts[:2]) if time_parts else 'less than 1ns'
            
            print(f'{message if message is not None else "Command"} '
                  f'{"aborted after" if aborted else "took"} '
                  f'{time_str}')
    else:
        yield  # no-op


def mergeDict(dict1, dict2):
    """ Merge dictionaries to collect info from several files"""
    dict3 = {**dict1, **dict2}
    for key, value in dict3.items():
        if key in dict1 and key in dict2:
                dict3[key] = [value , dict1[key]]
    return dict3

def flatten(lis):
     """ Recursively flattens a nested iterable"""
     for item in lis:
         if isinstance(item, Iterable) and not isinstance(item, str):
             for x in flatten(item):
                 yield x
         else:
             yield item

def pstdev(data):
    """Calculates the population standard deviation."""
    n = len(data)
    mean = sum(data)*1. / n  # mean
    var = sum(pow(x - mean, 2) for x in data) / n  # variance
    return math.sqrt(var)  # standard deviation

def find_polyA_motif(genome_seq, polyA_motif_list):
    """    
    Searches for the first occurrence of any polyA motif from a ranked list within a given genomic sequence.

    Args:
        genome_seq (str): The genomic sequence to search for polyA motifs. The sequence must already be oriented.
        polyA_motif_list (list of str): A ranked list of motifs to search for. The function will report the first motif found.

    Returns:
        tuple: A tuple containing:
            - polyA_motif (str): The first polyA motif found in the sequence. If no motif is found, returns 'NA'.
            - polyA_dist (int or str): The distance (in bases) upstream from the end of the sequence where the motif is found. If no motif is found, returns 'NA'.
            - found (str): 'TRUE' if a motif is found, otherwise 'FALSE'.
    """
    for motif in polyA_motif_list:
        i = genome_seq.find(motif)
        if i >= 0:
            return motif, -(len(genome_seq)-i-len(motif)+1), 'TRUE'
    return 'NA', 'NA', 'FALSE'


def get_files_from_dir(directory, extension):
    """ Get all files with a given extension from a directory or a file"""
    if os.path.isfile(directory):
        with open(directory) as f:
            return [line.strip() for line in f]  # Corrected strip method call
    else:
        return [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(extension)]
    

def find_closest_in_list(lst, pos):
    """
    Finds the closest value in a sorted list to a given position.

    Args:
        lst (list of int/float): A sorted list of numbers.
        pos (int/float): The position to find the closest value to.

    Returns:
        int/float: The difference between the closest value in the list and the given position.

    Example:
        >>> find_closest_in_list([1, 3, 5, 7], 4)
        -1
    """
    i = bisect.bisect_left(lst, pos)
    if i == 0:
        return lst[0]-pos
    elif i == len(lst):
        return lst[-1]-pos
    else:
        a, b = lst[i-1]-pos, lst[i]-pos
        if abs(a) < abs(b): return a
        else: return b
    
def alphanum_key(s):
    return [int(c) if c.isdigit() else c.lower() for c in re.split(r'(\d+)', s)]

    
def calculate_tss(strand, start0, end1):
    """
    Strand aware calculation of the middle of the peak in the bed file
    If the cage peak length is of 1 nucleotide, the average is not calculcated
    """
    if end1 - start0 > 1:
        tss0 = int((start0 + end1) / 2)
        if strand == '+':
            return tss0
        else:
            return tss0 + 1
    else:
        if strand == '+':
            return start0
        else:
            return end1
        
def read_gtf(file, attributes=["transcript_id"], keep_attributes=True):
    if keep_attributes:
        return pl.read_csv(file, separator="\t", comment_prefix="#", schema_overrides = {"seqname": pl.String}, has_header = False, new_columns=["seqname","source","feature","start","end","score","strand","frame","attributes"])\
            .with_columns(
                [pl.col("attributes").str.extract(rf'{attribute} "([^;]*)";').alias(attribute) for attribute in attributes]
                )
    else:
        return pl.read_csv(file, separator="\t", comment_prefix="#", schema_overrides = {"seqname": pl.String}, has_header = False, new_columns=["seqname","source","feature","start","end","score","strand","frame","attributes"])\
            .with_columns(
                [pl.col("attributes").str.extract(rf'{attribute} "([^;]*)";').alias(attribute) for attribute in attributes]
                ).drop("attributes")        