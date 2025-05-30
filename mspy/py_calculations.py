#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# -------------------------------------------------------------------------
# mspy/py_calculations.py
#
# Python re-implementation of functionalities from the original
# mspy/calculations.c module.
# -------------------------------------------------------------------------

import numpy as np
from scipy import signal, integrate, stats

# --- Data Structure Definitions (Conceptual) ---

class MBox:
    """
    Represents a bounding box.
    """
    def __init__(self, xmin, xmax, ymin, ymax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

    def __repr__(self):
        return f"MBox(xmin={self.xmin}, xmax={self.xmax}, ymin={self.ymin}, ymax={self.ymax})"

class MNoise:
    """
    Represents noise characteristics.
    """
    def __init__(self, level, width):
        self.level = level
        self.width = width

    def __repr__(self):
        return f"MNoise(level={self.level}, width={self.width})"

# --- Signal Printing ---

def array_print(arr):
    """
    Prints the contents of a 1D or 2D NumPy array.

    Args:
        arr (np.ndarray): The array to print.
    """
    print(arr)

# --- Basic Signal Operations ---

def signal_median(y_values_array):
    """
    Calculates the median of a NumPy array.

    Args:
        y_values_array (np.ndarray): Array of y-values.

    Returns:
        float: The median of the array. Returns np.nan for empty array.
    """
    if y_values_array.size == 0:
        return np.nan
    return np.median(y_values_array)

def signal_interpolate_y(x_known_array, y_known_array, x_target):
    """
    Performs linear interpolation for y at x_target.

    Args:
        x_known_array (np.ndarray): Array of known x-values.
        y_known_array (np.ndarray): Array of known y-values.
        x_target (float): The x-value at which to interpolate.

    Returns:
        float: The interpolated y-value.
    """
    return np.interp(x_target, x_known_array, y_known_array)

def signal_interpolate_x(x_known_array, y_known_array, y_target):
    """
    Performs linear interpolation for x at y_target.
    Assumes y_known_array is monotonic in the region of interest.

    Args:
        x_known_array (np.ndarray): Array of known x-values.
        y_known_array (np.ndarray): Array of known y-values.
        y_target (float): The y-value at which to interpolate.

    Returns:
        float: The interpolated x-value. Returns np.nan if y_target is outside y_known_array range
               or if y_known_array is not suitable for interpolation.
    """
    if y_known_array.size < 2:
        return np.nan
    # Ensure y_known_array is sorted to use np.interp correctly by swapping x and y
    try:
        # If y_known_array is monotonically increasing
        if np.all(np.diff(y_known_array) >= 0):
            return np.interp(y_target, y_known_array, x_known_array)
        # If y_known_array is monotonically decreasing
        elif np.all(np.diff(y_known_array) <= 0):
            return np.interp(y_target, y_known_array[::-1], x_known_array[::-1])
        else:
            # Non-monotonic array where interpolation is ambiguous
            return np.nan
    except Exception:
        return np.nan


# --- Signal Position and Size Functions ---

def signal_locate_x(x_array_sorted, x_value):
    """
    Finds the index where x_value would fit in a sorted x_array_sorted.

    Args:
        x_array_sorted (np.ndarray): A sorted array of x-values.
        x_value (float): The x-value to locate.

    Returns:
        int: The index where x_value would be inserted.
    """
    return np.searchsorted(x_array_sorted, x_value)

def signal_locate_max_y(y_array):
    """
    Finds the index of the maximum value in y_array.

    Args:
        y_array (np.ndarray): Array of y-values.

    Returns:
        int: The index of the maximum y-value. Returns -1 for an empty array.
    """
    if y_array.size == 0:
        return -1
    return np.argmax(y_array)

def signal_box(x_array, y_array):
    """
    Computes and returns the bounding box of the signal.

    Args:
        x_array (np.ndarray): Array of x-values.
        y_array (np.ndarray): Array of y-values.

    Returns:
        MBox: An MBox object representing the bounding box.
              Returns MBox with all NaNs if arrays are empty.
    """
    if x_array.size == 0 or y_array.size == 0:
        return MBox(np.nan, np.nan, np.nan, np.nan)
    return MBox(xmin=np.min(x_array), xmax=np.max(x_array),
                ymin=np.min(y_array), ymax=np.max(y_array))

# --- Peak-Related Functions ---

def signal_intensity(x_array, y_array, x_target):
    """
    Interpolates and returns the y-value (intensity) at x_target.

    Args:
        x_array (np.ndarray): Array of x-values.
        y_array (np.ndarray): Array of y-values.
        x_target (float): The x-value for which to find the intensity.

    Returns:
        float: The interpolated intensity at x_target.
    """
    return np.interp(x_target, x_array, y_array)

def signal_centroid(x_array, y_array, peak_x, peak_width, height_fraction=0.5):
    """
    Calculates the centroid of a peak in a signal.

    Args:
        x_array (np.ndarray): Array of x-values.
        y_array (np.ndarray): Array of y-values.
        peak_x (float): The x-coordinate of the peak maximum.
        peak_width (float): The approximate width of the peak.
        height_fraction (float, optional): Fraction of peak height to consider for centroid calculation.
                                           Defaults to 0.5 (FWHM).

    Returns:
        float: The x-coordinate of the peak centroid. Returns np.nan if peak cannot be processed.
    """
    if x_array.size == 0 or y_array.size == 0:
        return np.nan

    peak_y = np.interp(peak_x, x_array, y_array)
    threshold_y = peak_y * height_fraction

    # Define a window around the peak_x to look for centroid
    window_half_width = peak_width * 1.5 # A bit wider than the peak_width for safety
    min_x_window = peak_x - window_half_width
    max_x_window = peak_x + window_half_width

    indices_in_window = np.where((x_array >= min_x_window) & (x_array <= max_x_window) & (y_array >= threshold_y))[0]

    if indices_in_window.size == 0:
        return np.nan # No points above threshold in window

    selected_x = x_array[indices_in_window]
    selected_y = y_array[indices_in_window]

    # Centroid calculation: sum(x_i * y_i) / sum(y_i)
    # Subtracting the threshold makes it more robust to baseline issues
    weighted_sum_x = np.sum(selected_x * (selected_y - threshold_y))
    sum_weights = np.sum(selected_y - threshold_y)

    if sum_weights == 0:
        return np.nan # Avoid division by zero

    return weighted_sum_x / sum_weights

def signal_width(x_array, y_array, peak_x, peak_y, height_fraction=0.5):
    """
    Computes the width of the peak at a specified height_fraction of peak_y.

    Args:
        x_array (np.ndarray): Array of x-values.
        y_array (np.ndarray): Array of y-values.
        peak_x (float): The x-coordinate of the peak maximum.
        peak_y (float): The intensity of the peak maximum.
        height_fraction (float, optional): Fraction of peak_y at which to measure width. Defaults to 0.5 (FWHM).

    Returns:
        float: The width of the peak. Returns np.nan if width cannot be determined.
    """
    if x_array.size < 2 or y_array.size < 2:
        return np.nan

    target_height = peak_y * height_fraction
    
    # Find indices where y_array crosses target_height
    # Consider points to the left and right of the peak maximum separately
    peak_idx = np.searchsorted(x_array, peak_x)
    if peak_idx >= x_array.size: # peak_x is beyond the rightmost x_array point
        peak_idx = x_array.size -1
    elif x_array[peak_idx] != peak_x and peak_idx > 0: # Adjust if peak_x is not exactly in x_array
        if abs(x_array[peak_idx] - peak_x) > abs(x_array[peak_idx-1] - peak_x):
            peak_idx -=1
            
    # Left side
    x_left = np.nan
    left_indices = np.where(y_array[:peak_idx+1] >= target_height)[0]
    if left_indices.size > 0:
        # Last point above or at target height on the left
        idx1_left = left_indices[-1]
        if idx1_left > 0: # We need a point before it to interpolate
            x1, y1 = x_array[idx1_left-1], y_array[idx1_left-1]
            x2, y2 = x_array[idx1_left], y_array[idx1_left]
            if y1 < target_height < y2 : # Crosses from below
                 x_left = x1 + (x2 - x1) * (target_height - y1) / (y2 - y1)
            elif y2 == target_height:
                 x_left = x2
            elif idx1_left +1 < peak_idx+1 and y_array[idx1_left+1] < target_height and y_array[idx1_left] >= target_height : # Crosses from above (flat top)
                 x_left = x_array[idx1_left]


    # Right side
    x_right = np.nan
    right_indices = np.where(y_array[peak_idx:] >= target_height)[0]
    if right_indices.size > 0:
        # First point above or at target height on the right (relative to peak_idx)
        idx1_right = right_indices[0] + peak_idx
        if idx1_right < y_array.size -1 : # We need a point after it to interpolate
            x1, y1 = x_array[idx1_right], y_array[idx1_right]
            x2, y2 = x_array[idx1_right+1], y_array[idx1_right+1]
            if y1 > target_height > y2 : # Crosses from above
                x_right = x1 + (x2 - x1) * (target_height - y1) / (y2 - y1)
            elif y1 == target_height:
                x_right = x1
            elif idx1_right > 0 and y_array[idx1_right-1] < target_height and y_array[idx1_right] >= target_height: # Crosses from below (flat top)
                 x_right = x_array[idx1_right]


    if not np.isnan(x_left) and not np.isnan(x_right):
        return x_right - x_left
    else:
        return np.nan

def signal_area(x_array, y_array):
    """
    Calculates the area under the signal curve using the trapezoidal rule.

    Args:
        x_array (np.ndarray): Array of x-values.
        y_array (np.ndarray): Array of y-values.

    Returns:
        float: The area under the curve. Returns 0 for empty or single-point arrays.
    """
    if x_array.size < 2 or y_array.size < 2:
        return 0.0
    return integrate.trapz(y_array, x_array)

def signal_noise(y_array):
    """
    Estimates noise level and width using Median Absolute Deviation (MAD).

    Args:
        y_array (np.ndarray): Array of y-values.

    Returns:
        MNoise: An MNoise object with 'level' and 'width'.
                Level is MAD, width is also MAD (can be adapted if a different width metric is needed).
                Returns MNoise(np.nan, np.nan) for empty array.
    """
    if y_array.size == 0:
        return MNoise(np.nan, np.nan)
    
    # SciPy's median_abs_deviation by default scales by 1/stats.norm.ppf(0.75) which is approx 1.4826
    # to estimate the standard deviation of a normally distributed variable.
    # The original C code used a factor of 1.0, so we set scale='normal' to get the factor,
    # or we can calculate MAD manually and multiply by 1.4826 if we want to match that specific scaling.
    # For simplicity and robustness, using scipy.stats.median_abs_deviation directly is good.
    # The 'width' was not clearly defined in the C context, using MAD level as a proxy for now.
    # If a different "width" of noise is needed, this part might need adjustment.
    
    mad_level = stats.median_abs_deviation(y_array, scale='normal') # Consistent with 1.4826*median(abs(dev))
    # If the original C code's 'width' had a different meaning, this needs adjustment.
    # For now, using the same MAD value for width as a placeholder.
    return MNoise(level=mad_level, width=mad_level)


def signal_local_maxima(y_array, **kwargs):
    """
    Identifies and returns local maxima (peaks) in y_array.
    Uses scipy.signal.find_peaks.

    Args:
        y_array (np.ndarray): Array of y-values.
        **kwargs: Additional keyword arguments to pass to scipy.signal.find_peaks
                  (e.g., height, prominence, width, distance).

    Returns:
        tuple: A tuple containing:
            - peak_indices (np.ndarray): Indices of the peaks in y_array.
            - properties (dict): Dictionary of properties associated with each peak.
    """
    if y_array.size == 0:
        return np.array([], dtype=int), {}
    peak_indices, properties = signal.find_peaks(y_array, **kwargs)
    return peak_indices, properties

# --- Signal Manipulation Functions ---

def signal_crop(x_array, y_array, x_min_crop, x_max_crop):
    """
    Extracts a segment of the signal between x_min_crop and x_max_crop.
    Interpolates start/end points if crop boundaries fall between original points.

    Args:
        x_array (np.ndarray): Array of x-values.
        y_array (np.ndarray): Array of y-values.
        x_min_crop (float): Minimum x-value for cropping.
        x_max_crop (float): Maximum x-value for cropping.

    Returns:
        tuple: (x_cropped, y_cropped) as NumPy arrays.
               Returns empty arrays if crop range is invalid or no data falls within it.
    """
    if x_array.size == 0 or y_array.size == 0 or x_min_crop >= x_max_crop:
        return np.array([]), np.array([])

    # Find indices for cropping range
    start_idx = np.searchsorted(x_array, x_min_crop, side='left')
    end_idx = np.searchsorted(x_array, x_max_crop, side='right')

    # Create lists to build the new cropped arrays
    x_cropped_list = []
    y_cropped_list = []

    # Interpolate and add the start point if x_min_crop is within the original range
    # and not an exact match to an existing x-value (unless it's the very first point)
    if x_min_crop > x_array[0] and x_min_crop < x_array[-1]:
        y_start_interp = np.interp(x_min_crop, x_array, y_array)
        x_cropped_list.append(x_min_crop)
        y_cropped_list.append(y_start_interp)

    # Add original points within the crop range
    # Adjust start_idx if x_min_crop is to the left of the first point in the selected slice
    if start_idx < x_array.size and x_array[start_idx] < x_min_crop :
        actual_start_idx = start_idx +1
    else:
        actual_start_idx = start_idx
        
    # Adjust end_idx similarly
    if end_idx > 0 and x_array[end_idx-1] > x_max_crop:
        actual_end_idx = end_idx -1
    else:
        actual_end_idx = end_idx


    x_cropped_list.extend(x_array[actual_start_idx:actual_end_idx])
    y_cropped_list.extend(y_array[actual_start_idx:actual_end_idx])
    
    # Interpolate and add the end point if x_max_crop is within the original range
    # and not an exact match to an existing x-value (unless it's the very last point)
    if x_max_crop > x_array[0] and x_max_crop < x_array[-1]:
        y_end_interp = np.interp(x_max_crop, x_array, y_array)
        x_cropped_list.append(x_max_crop)
        y_cropped_list.append(y_end_interp)
        
    # Ensure sorted x values and remove duplicates that might arise from interpolation
    if x_cropped_list:
        x_cropped = np.array(x_cropped_list)
        y_cropped = np.array(y_cropped_list)
        sorted_indices = np.argsort(x_cropped)
        x_cropped = x_cropped[sorted_indices]
        y_cropped = y_cropped[sorted_indices]
        
        # Remove duplicates keeping the first occurrence (important if x_min/max_crop coincided with existing points)
        unique_x, unique_indices = np.unique(x_cropped, return_index=True)
        if unique_x.size < x_cropped.size: # Only filter if duplicates exist
            x_cropped = x_cropped[unique_indices]
            y_cropped = y_cropped[unique_indices]
        return x_cropped, y_cropped
    else: # If no points were added (e.g., crop range outside data)
        return np.array([]), np.array([])


def signal_offset(x_array, y_array, x_offset, y_offset):
    """
    Offsets the signal by adding x_offset to x-values and y_offset to y-values.

    Args:
        x_array (np.ndarray): Array of x-values.
        y_array (np.ndarray): Array of y-values.
        x_offset (float): Value to add to x-values.
        y_offset (float): Value to add to y-values.

    Returns:
        tuple: (x_offsetted, y_offsetted) as NumPy arrays.
    """
    return x_array + x_offset, y_array + y_offset

def signal_multiply(x_array, y_array, x_factor, y_factor):
    """
    Multiplies the signal by x_factor for x-values and y_factor for y-values.

    Args:
        x_array (np.ndarray): Array of x-values.
        y_array (np.ndarray): Array of y-values.
        x_factor (float): Factor to multiply x-values by.
        y_factor (float): Factor to multiply y-values by.

    Returns:
        tuple: (x_multiplied, y_multiplied) as NumPy arrays.
    """
    return x_array * x_factor, y_array * y_factor

def signal_normalize(y_array):
    """
    Normalizes y-values so the maximum becomes 1.0.
    If max is 0 or array is empty, returns a copy of the original y_array (or empty).

    Args:
        y_array (np.ndarray): Array of y-values.

    Returns:
        np.ndarray: Normalized y_array.
    """
    if y_array.size == 0:
        return np.array([])
    max_y = np.max(y_array)
    if max_y == 0:
        return y_array.copy()  # Avoid division by zero, return as is or zeros
    return y_array / max_y

def signal_smooth_ma(y_array, window_size, cycles=1):
    """
    Smooths y_array using a moving average filter over multiple cycles.

    Args:
        y_array (np.ndarray): Array of y-values to smooth.
        window_size (int): Size of the moving average window. Must be an odd integer.
        cycles (int, optional): Number of smoothing cycles. Defaults to 1.

    Returns:
        np.ndarray: The smoothed y_array.
    """
    if window_size <= 0:
        raise ValueError("window_size must be a positive integer.")
    if window_size % 2 == 0:
        # Adjust to make it odd, or raise error. For simplicity, let's adjust.
        # print("Warning: window_size should be odd for 'same' mode without phase shift. Adjusting.")
        # window_size +=1 # This might not be desired, better to enforce odd or handle padding carefully.
        # For 'same' mode, an even window will cause a half-sample shift if not handled.
        # However, numpy.convolve with 'same' handles this by making the output length same as input.
        pass


    smoothed_y = y_array.copy()
    if window_size == 1: # No smoothing if window is 1
        return smoothed_y

    for _ in range(cycles):
        if smoothed_y.size < window_size: # Cannot convolve if array is smaller than window
            # Return array as is or handle as an error/warning
            # print(f"Warning: array size ({smoothed_y.size}) is smaller than window_size ({window_size}). Skipping smoothing.")
            return smoothed_y
        window = np.ones(window_size) / float(window_size)
        smoothed_y = np.convolve(smoothed_y, window, mode='same')
    return smoothed_y

# --- Formula Composition (Placeholder, if needed from calculations.c) ---
# The calculations.c file had a 'formula_composition' function, which seems to be
# related to generating molecular formulas based on mass. This is a complex
# combinatorial problem. A direct Python port of such a C function would be
# non-trivial and likely slow if not carefully optimized or if it relies on
# specific C optimizations.
#
# For now, we'll assume that the primary need from calculations.c was for the
# signal processing functions. If formula_composition is indeed required,
# it would need a dedicated implementation, possibly leveraging libraries
# like a Python-based isotopic distribution calculator or a more specialized
# chemical formula generator.
#
# A placeholder function or a note indicating its absence:

def formula_composition(minimum_counts, maximum_counts, element_masses,
                        low_mass_target, high_mass_target, max_results_limit):
    """
    Placeholder for a function that generates molecular compositions.
    This functionality was present in the C extension but requires a more
    complex implementation in Python, possibly using specialized libraries.

    Args:
        minimum_counts (tuple): Minimum atom counts for each element.
        maximum_counts (tuple): Maximum atom counts for each element.
        element_masses (tuple): Masses of the elements (sorted reverse).
        low_mass_target (float): Low mass limit for generated formulas.
        high_mass_target (float): High mass limit for generated formulas.
        max_results_limit (int): Maximum number of formula combinations.

    Returns:
        list: A list of compositions (e.g., lists of counts).
              Currently returns an empty list as it's a placeholder.
    """
    print("Warning: formula_composition is a placeholder and not fully implemented in py_calculations.py.")
    # This would be a complex combinatorial search.
    # Example of how one might start thinking about it (very simplified):
    # results = []
    # from itertools import product
    # count_ranges = [range(min_c, max_c + 1) for min_c, max_c in zip(minimum_counts, maximum_counts)]
    # for combo in product(*count_ranges):
    #     current_mass = sum(c * m for c, m in zip(combo, element_masses))
    #     if low_mass_target <= current_mass <= high_mass_target:
    #         results.append(list(combo))
    #         if len(results) >= max_results_limit:
    #             break
    # return results
    # return [] # Original placeholder

    results = []
    num_elements = len(masses)
    
    # current_composition will be mutated by the recursive helper
    current_composition_counts = [0] * num_elements

    # Pre-calculate min and max possible masses for remaining elements for pruning
    min_mass_remaining_at_idx = [0.0] * (num_elements + 1)
    max_mass_remaining_at_idx = [0.0] * (num_elements + 1)
    for i in range(num_elements - 1, -1, -1):
        min_mass_remaining_at_idx[i] = minimum_counts[i] * masses[i] + min_mass_remaining_at_idx[i+1]
        max_mass_remaining_at_idx[i] = maximum_counts[i] * masses[i] + max_mass_remaining_at_idx[i+1]


    def generate_recursive(element_idx, current_mass_sum):
        if len(results) >= max_results_limit:
            return

        if element_idx == num_elements: # Base case: all elements have been assigned a count
            if lo_mass_target <= current_mass_sum <= high_mass_target:
                results.append(list(current_composition_counts))
            return

        # Iterate for current element element_idx
        original_count_for_elem = current_composition_counts[element_idx] # Should be 0 if coming from parent level correctly

        for count in range(minimum_counts[element_idx], maximum_counts[element_idx] + 1):
            current_composition_counts[element_idx] = count
            new_mass_intermediate = current_mass_sum + count * masses[element_idx]

            # Pruning condition 1: If current mass already exceeds hi_mass
            if new_mass_intermediate > high_mass_target and masses[element_idx] > 0: # check for positive mass element
                 # Since masses are sorted descending, further counts for this element or deeper elements will also exceed.
                 # However, this only works if all masses[i] >= 0. If negative masses are allowed (unlikely for elements), this logic changes.
                 # Assuming positive masses.
                current_composition_counts[element_idx] = original_count_for_elem # backtrack for this element
                return # Prune this path and parent counts for this element_idx

            # Pruning condition 2: If current_mass + min_possible_for_remaining > hi_mass
            # Check if new_mass_intermediate plus the minimum possible mass from *subsequent* elements would exceed hi_mass
            if element_idx + 1 < num_elements:
                min_mass_from_next_elements = min_mass_remaining_at_idx[element_idx+1]
                if new_mass_intermediate + min_mass_from_next_elements > high_mass_target:
                    # If even by picking the smallest counts for all subsequent elements, we exceed the target,
                    # then increasing the count of the current element (element_idx) will only make it worse.
                    # So, we can stop trying further counts for element_idx.
                    current_composition_counts[element_idx] = original_count_for_elem # backtrack
                    return # Prune this path and parent counts for this element_idx
            elif new_mass_intermediate > high_mass_target: # Last element, check directly
                    current_composition_counts[element_idx] = original_count_for_elem # backtrack
                    return


            # Pruning condition 3: If current_mass + max_possible_for_remaining < lo_mass
            # Check if new_mass_intermediate plus the maximum possible mass from *subsequent* elements is still less than lo_mass
            if element_idx + 1 < num_elements:
                max_mass_from_next_elements = max_mass_remaining_at_idx[element_idx+1]
                if new_mass_intermediate + max_mass_from_next_elements < lo_mass_target:
                    # If even by picking the largest counts for all subsequent elements, we are still below the target,
                    # then the current count for element_idx is too small. So, continue to the next count for element_idx.
                    # No need to backtrack current_composition_counts[element_idx] as it's reset by the loop.
                    continue 
            elif new_mass_intermediate < lo_mass_target: # Last element, check directly
                    continue


            # If not pruned, go to the next element
            generate_recursive(element_idx + 1, new_mass_intermediate)
            
            if len(results) >= max_results_limit:
                current_composition_counts[element_idx] = original_count_for_elem # backtrack
                return
        
        current_composition_counts[element_idx] = original_count_for_elem # Backtrack after trying all counts for this element

    generate_recursive(0, 0.0)
    return results

# --- Peak Shape Models (Gaussian, Lorentzian etc. - if needed from calculations.c) ---
# The C file also had functions for generating peak shapes (gaussian, lorentzian).
# These are often used for creating synthetic spectra or for fitting.
# For basic profile generation, these might not be directly called if using
# higher-level functions from SciPy or Matplotlib for plotting/simulation.
# However, if direct generation of a single peak shape is needed:

def signal_gaussian(x_center, y_min, y_max, fwhm, num_points=500):
    """
    Generates a Gaussian peak shape.

    Args:
        x_center (float): X-coordinate of the peak center.
        y_min (float): Minimum y-value (baseline).
        y_max (float): Maximum y-value (peak height above baseline).
        fwhm (float): Full Width at Half Maximum of the peak.
        num_points (int, optional): Number of points to generate for the peak. Defaults to 500.

    Returns:
        tuple: (x_coords, y_coords) for the Gaussian peak.
    """
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    # Generate x points covering roughly +/- 3-4 sigmas from the center
    # A common choice is to cover where the peak is significant
    x_range = np.linspace(x_center - 3.5 * sigma, x_center + 3.5 * sigma, num_points)
    y_coords = y_min + (y_max - y_min) * np.exp(-((x_range - x_center)**2) / (2 * sigma**2))
    return x_range, y_coords

def signal_lorentzian(x_center, y_min, y_max, fwhm, num_points=500):
    """
    Generates a Lorentzian peak shape.

    Args:
        x_center (float): X-coordinate of the peak center.
        y_min (float): Minimum y-value (baseline).
        y_max (float): Maximum y-value (peak height above baseline).
        fwhm (float): Full Width at Half Maximum of the peak.
        num_points (int, optional): Number of points to generate for the peak. Defaults to 500.

    Returns:
        tuple: (x_coords, y_coords) for the Lorentzian peak.
    """
    gamma = fwhm / 2.0 # gamma is HWHM for Lorentzian
    # Generate x points covering a wider range for Lorentzian due to heavier tails
    x_range = np.linspace(x_center - 10 * gamma, x_center + 10 * gamma, num_points)
    y_coords = y_min + (y_max - y_min) * (gamma**2 / ((x_range - x_center)**2 + gamma**2))
    return x_range, y_coords

def signal_gausslorentzian(x_center, y_min, y_max, fwhm, num_points=500, gaussian_fraction=0.5):
    """
    Generates a pseudo-Voigt peak shape (linear combination of Gaussian and Lorentzian).

    Args:
        x_center (float): X-coordinate of the peak center.
        y_min (float): Minimum y-value (baseline).
        y_max (float): Maximum y-value (peak height above baseline).
        fwhm (float): Full Width at Half Maximum of the peak.
        num_points (int, optional): Number of points to generate for the peak. Defaults to 500.
        gaussian_fraction (float, optional): Fraction of Gaussian component (0 to 1). Defaults to 0.5.

    Returns:
        tuple: (x_coords, y_coords) for the pseudo-Voigt peak.
    """
    _, gauss_y = signal_gaussian(x_center, 0, 1, fwhm, num_points) # Normalized height
    x_range_lor, lor_y = signal_lorentzian(x_center, 0, 1, fwhm, num_points) # Normalized height

    # Assuming x_range from lorentzian is suitable (often wider)
    # Or, one could define a common x_range for both. For simplicity, using Lorentzian's.
    # Re-evaluate gaussian on the lorentzian's x_range if they differ significantly.
    # For this implementation, let's assume the x_range from lorentzian is adequate.
    # If num_points is the same, and fwhm leads to similar effective widths, this is often okay.
    
    # A more robust way is to define a common x_range first:
    sigma = fwhm / (2 * np.sqrt(2 * np.log(2)))
    gamma = fwhm / 2.0
    # Use a range that covers both well, e.g., based on the broader Lorentzian tails
    x_range_common = np.linspace(x_center - 10 * gamma, x_center + 10 * gamma, num_points)
    
    gauss_y_common = np.exp(-((x_range_common - x_center)**2) / (2 * sigma**2))
    lor_y_common = (gamma**2 / ((x_range_common - x_center)**2 + gamma**2))

    y_coords = y_min + (y_max - y_min) * (gaussian_fraction * gauss_y_common + (1 - gaussian_fraction) * lor_y_common)
    # Normalization might be needed if the sum of fractions isn't strictly 1 or if peaks aren't normalized to 1 before combining
    # The above assumes that both gauss_y_common and lor_y_common peak at 1 at x_center
    
    return x_range_common, y_coords

# Example of how the C extension's profile generation might be used (conceptual)
# The actual C code uses calculations.signal_profile and calculations.signal_profile_to_raster
# These would take a list of peaks (mass, intensity, fwhm) and generate a full profile.
# This is more complex than generating single peaks.
# A Python equivalent would involve summing up individual peak shapes.

def signal_profile(peaks_data, points_per_fwhm=10, noise_level=0, model='gaussian'):
    """
    Generates a profile spectrum from a list of peaks.

    Args:
        peaks_data (list of lists/tuples): Each item: [mz, intensity, fwhm].
        points_per_fwhm (int): Number of points to generate across FWHM of each peak.
        noise_level (float): Amplitude of random noise to add.
        model (str): 'gaussian', 'lorentzian', or 'gausslorentzian'.

    Returns:
        np.ndarray: A 2D array where rows are [x, y] coordinates of the profile.
                    Returns an empty array if peaks_data is empty.
    """
    if not peaks_data:
        return np.array([]).reshape(0, 2)

    all_x = []
    all_y = []

    for mz, intensity, fwhm_val in peaks_data:
        if fwhm_val <= 0: continue # Skip invalid peaks

        num_points_peak = int(points_per_fwhm * 5) # Heuristic: cover ~2.5 FWHM on each side
        if num_points_peak < 5: num_points_peak = 5


        if model == 'gaussian':
            x_peak, y_peak = signal_gaussian(mz, 0, intensity, fwhm_val, num_points_peak)
        elif model == 'lorentzian':
            x_peak, y_peak = signal_lorentzian(mz, 0, intensity, fwhm_val, num_points_peak)
        elif model == 'gausslorentzian': # Assuming 0.5 gaussian fraction by default
            x_peak, y_peak = signal_gausslorentzian(mz, 0, intensity, fwhm_val, num_points_peak)
        else:
            raise ValueError(f"Unknown peak model: {model}")
        
        all_x.extend(x_peak)
        all_y.extend(y_peak)

    if not all_x:
        return np.array([]).reshape(0,2)

    # Sort by x values and combine overlapping points by summing y values
    # This is a simplified way; proper merging might need resampling onto a common grid
    sorted_indices = np.argsort(all_x)
    x_sorted = np.array(all_x)[sorted_indices]
    y_sorted = np.array(all_y)[sorted_indices]

    # For simplicity, we'll use a common grid based on min/max x and a resolution
    # A more advanced approach would use np.histogram or similar for binning and summing
    if x_sorted.size == 0:
      return np.array([]).reshape(0,2)
      
    # Determine overall range and a suitable number of points for the final profile
    # This resolution needs to be fine enough to capture individual peak shapes
    # A simple approach: use a multiple of the total number of points generated for all peaks,
    # or determine based on smallest FWHM.
    # For now, let's define a reasonable number of points for the output profile.
    # This part is tricky to match C behavior without knowing its exact raster generation.
    
    # Let's define a global raster for the output for simplicity,
    # though the C code might do adaptive rastering.
    # This is a common approach in Python plotting libraries.
    
    min_x_overall = np.min(x_sorted)
    max_x_overall = np.max(x_sorted)
    
    # Heuristic for number of points in the final profile:
    # Try to get at least 'points_per_fwhm' for the narrowest peak.
    min_fwhm = min(p[2] for p in peaks_data if p[2] > 0) if peaks_data else 0.1
    if min_fwhm <= 0: min_fwhm = 0.001 # Avoid division by zero for very narrow/zero fwhm
    
    num_final_points = int( (max_x_overall - min_x_overall) / (min_fwhm / points_per_fwhm) )
    if num_final_points < 2 * len(peaks_data): # Ensure at least a few points per peak on average
        num_final_points = 2 * len(peaks_data)
    if num_final_points < 100: num_final_points = 100 # Minimum number of points
    if num_final_points > 20000: num_final_points = 20000 # Cap max points for performance

    profile_x_coords = np.linspace(min_x_overall, max_x_overall, num_final_points)
    profile_y_coords = np.zeros_like(profile_x_coords)

    # Add each peak to the common raster
    for mz, intensity, fwhm_val in peaks_data:
        if fwhm_val <= 0: continue

        if model == 'gaussian':
            sigma = fwhm_val / (2 * np.sqrt(2 * np.log(2)))
            y_add = intensity * np.exp(-((profile_x_coords - mz)**2) / (2 * sigma**2))
        elif model == 'lorentzian':
            gamma = fwhm_val / 2.0
            y_add = intensity * (gamma**2 / ((profile_x_coords - mz)**2 + gamma**2))
        elif model == 'gausslorentzian':
            sigma = fwhm_val / (2 * np.sqrt(2 * np.log(2)))
            gamma = fwhm_val / 2.0
            gauss_comp = 0.5 * intensity * np.exp(-((profile_x_coords - mz)**2) / (2 * sigma**2))
            lor_comp = 0.5 * intensity * (gamma**2 / ((profile_x_coords - mz)**2 + gamma**2))
            y_add = gauss_comp + lor_comp
        else: # Should not happen given earlier check
            y_add = np.zeros_like(profile_x_coords)
            
        profile_y_coords += y_add

    if noise_level > 0:
        profile_y_coords += np.random.normal(0, noise_level * np.max(profile_y_coords) if np.max(profile_y_coords) > 0 else noise_level, size=profile_y_coords.shape)
        profile_y_coords[profile_y_coords < 0] = 0 # Ensure noise doesn't make intensity negative

    return np.vstack((profile_x_coords, profile_y_coords)).T


def signal_profile_to_raster(peaks_data, raster_x, noise_level=0, model_shape=0):
    """
    Generates a profile spectrum from a list of peaks onto a predefined x-axis raster.

    Args:
        peaks_data (np.ndarray): 2D array where rows are [mz, intensity, fwhm].
        raster_x (np.ndarray): 1D array of x-values for the output raster.
        noise_level (float): Amplitude of random noise to add (relative to max intensity if >0, else absolute).
        model_shape (int): 0 for Gaussian, 1 for Lorentzian, 2 for Gauss-Lorentzian (Pseudo-Voigt).

    Returns:
        np.ndarray: A 2D array where rows are [x, y] coordinates of the profile on the raster.
    """
    if peaks_data.ndim != 2 or peaks_data.shape[1] != 3:
        raise ValueError("peaks_data must be a 2D array with columns [mz, intensity, fwhm]")
    if raster_x.ndim != 1:
        raise ValueError("raster_x must be a 1D array.")
    if raster_x.size == 0:
        return np.array([]).reshape(0,2)

    profile_y_coords = np.zeros_like(raster_x, dtype=float)

    for mz, intensity, fwhm_val in peaks_data:
        if fwhm_val <= 0 or intensity <=0: continue

        if model_shape == 0: # Gaussian
            sigma = fwhm_val / (2 * np.sqrt(2 * np.log(2)))
            if sigma == 0: continue
            y_add = intensity * np.exp(-((raster_x - mz)**2) / (2 * sigma**2))
        elif model_shape == 1: # Lorentzian
            gamma = fwhm_val / 2.0
            if gamma == 0: continue
            y_add = intensity * (gamma**2 / ((raster_x - mz)**2 + gamma**2))
        elif model_shape == 2: # Gauss-Lorentzian (Pseudo-Voigt, assuming 0.5 fraction)
            sigma = fwhm_val / (2 * np.sqrt(2 * np.log(2)))
            gamma = fwhm_val / 2.0
            if sigma == 0 or gamma == 0: continue # Avoid division by zero if FWHM is tiny
            gauss_comp = 0.5 * intensity * np.exp(-((raster_x - mz)**2) / (2 * sigma**2))
            lor_comp = 0.5 * intensity * (gamma**2 / ((raster_x - mz)**2 + gamma**2))
            y_add = gauss_comp + lor_comp
        else:
            raise ValueError(f"Unknown peak model_shape: {model_shape}")
            
        profile_y_coords += y_add
    
    if noise_level > 0:
        max_intensity = np.max(profile_y_coords) if profile_y_coords.size > 0 and np.max(profile_y_coords) > 0 else 1.0
        # Add noise relative to max signal or absolute if max signal is 0
        actual_noise_std = noise_level * max_intensity if max_intensity > 0 else noise_level
        profile_y_coords += np.random.normal(0, actual_noise_std, size=profile_y_coords.shape)
        profile_y_coords[profile_y_coords < 0] = 0 # Ensure noise doesn't make intensity negative

    return np.vstack((raster_x, profile_y_coords)).T
