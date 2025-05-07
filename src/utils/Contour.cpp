/**
 * Contour.cpp
 * 
 * Implementation of the Contour class for efficient packing in the ASF-B*-tree
 * placement algorithm. The contour data structure represents the skyline profile
 * of the currently placed modules, allowing for O(log n) height queries and updates.
 */

#include "Contour.hpp"
#include <algorithm>
#include <cassert>

/**
 * Constructor
 */
Contour::Contour() {
    // Initialize with empty contour
}

/**
 * Copy constructor
 */
Contour::Contour(const Contour& other) : contourMap(other.contourMap) {
    // Copy the contour map from the other contour
}

/**
 * Destructor
 */
Contour::~Contour() {
    // No manual cleanup needed
}

/**
 * Clears the contour
 */
void Contour::clear() {
    contourMap.clear();
}

/**
 * Adds a segment to the contour
 * 
 * This function updates the contour to incorporate a new module or segment.
 * It ensures the contour remains a valid skyline by merging overlapping segments.
 * 
 * Time complexity: O(log n) for insertion + O(k) for merging where k is the
 * number of contour segments affected.
 */
void Contour::addSegment(int start, int end, int height) {
    if (start >= end) return;  // Invalid segment
    
    // Find the first contour segment whose start point is >= our start
    auto it = contourMap.lower_bound(start);
    
    // Check if we need to update the previous segment
    if (it != contourMap.begin()) {
        auto prev = std::prev(it);
        if (prev->second > height) {
            // If previous segment is higher, we need to insert our start
            contourMap[start] = height;
        } else if (prev->second < height) {
            // If previous segment is lower, do nothing (will be handled below)
        } else {
            // If same height, we can extend the previous segment
            start = prev->first;
            contourMap.erase(prev);
        }
    }
    
    // Remove all segments that are completely covered by our new segment
    while (it != contourMap.end() && it->first < end) {
        auto next = std::next(it);
        
        if (it->second <= height) {
            // If current segment is not higher than our new segment, remove it
            contourMap.erase(it);
        } else {
            // If current segment is higher, we need to keep it
            if (it->first == start) {
                // If it starts at our start, do nothing
            } else {
                // Otherwise, we need to add a new point at our start
                contourMap[start] = height;
            }
            
            // Skip to the end of this higher segment
            break;
        }
        
        it = next;
    }
    
    // Add our end point
    if (it == contourMap.end() || it->first > end) {
        contourMap[end] = it != contourMap.end() ? it->second : 0;
    }
    
    // Finally, add our start point
    contourMap[start] = height;
}

/**
 * Gets the height of the contour at a specific range
 * 
 * This function returns the maximum height of the contour within
 * the specified coordinate range.
 * 
 * Time complexity: O(log n + k) where k is the number of segments in the range
 */
int Contour::getHeight(int start, int end) const {
    if (start >= end) return 0;  // Invalid range
    if (contourMap.empty()) return 0;  // Empty contour
    
    // Find the first segment that starts at or before our query start
    auto it = contourMap.upper_bound(start);
    if (it != contourMap.begin()) {
        --it;  // Move to the segment that contains our start
    }
    
    int maxHeight = 0;
    
    // Iterate through all segments that overlap with our range
    while (it != contourMap.end() && it->first < end) {
        maxHeight = std::max(maxHeight, it->second);
        ++it;
    }
    
    return maxHeight;
}

/**
 * Gets all contour segments
 * 
 * This function converts the internal representation to a vector
 * of ContourSegment objects for easier processing.
 * 
 * Time complexity: O(n) where n is the number of segments
 */
std::vector<ContourSegment> Contour::getSegments() const {
    std::vector<ContourSegment> segments;
    
    if (contourMap.empty()) {
        return segments;
    }
    
    // Convert map representation to segments
    auto it = contourMap.begin();
    auto prev = it;
    ++it;
    
    while (it != contourMap.end()) {
        segments.emplace_back(prev->first, it->first, prev->second);
        prev = it;
        ++it;
    }
    
    return segments;
}

/**
 * Merges this contour with another contour
 * 
 * This function combines two contours by taking the maximum height
 * at each coordinate.
 * 
 * Time complexity: O(n + m) where n and m are the sizes of the two contours
 */
void Contour::merge(const Contour& other) {
    // Get all the breakpoints from both contours
    std::vector<int> breakpoints;
    
    for (const auto& point : contourMap) {
        breakpoints.push_back(point.first);
    }
    
    for (const auto& point : other.contourMap) {
        breakpoints.push_back(point.first);
    }
    
    // Sort and remove duplicates
    std::sort(breakpoints.begin(), breakpoints.end());
    breakpoints.erase(std::unique(breakpoints.begin(), breakpoints.end()), breakpoints.end());
    
    // Create a new contour map
    std::map<int, int> newContourMap;
    
    // For each breakpoint, take the maximum height from both contours
    for (size_t i = 0; i < breakpoints.size(); ++i) {
        int coordinate = breakpoints[i];
        int height1 = getHeight(coordinate, coordinate + 1);
        int height2 = other.getHeight(coordinate, coordinate + 1);
        
        newContourMap[coordinate] = std::max(height1, height2);
    }
    
    // Update our contour map
    contourMap = std::move(newContourMap);
}

/**
 * Gets the maximum coordinate value in the contour
 */
int Contour::getMaxCoordinate() const {
    if (contourMap.empty()) {
        return 0;
    }
    
    return contourMap.rbegin()->first;
}

/**
 * Gets the maximum height value in the contour
 */
int Contour::getMaxHeight() const {
    int maxHeight = 0;
    
    for (const auto& point : contourMap) {
        maxHeight = std::max(maxHeight, point.second);
    }
    
    return maxHeight;
}

/**
 * Checks if the contour is empty
 */
bool Contour::isEmpty() const {
    return contourMap.empty();
}

/**
 * Gets the contour map for direct access
 */
const std::map<int, int>& Contour::getContourMap() const {
    return contourMap;
}