/**
 * Contour.hpp
 * 
 * This file defines the Contour class used for efficient packing in the ASF-B*-tree
 * placement algorithm. The contour data structure represents the skyline profile
 * of the currently placed modules, allowing for O(log n) height queries and updates.
 */

#pragma once

#include <map>
#include <vector>
#include <utility>
#include <limits>
#include <memory>

/**
 * Contour segment representing a horizontal or vertical segment in the placement
 */
struct ContourSegment {
    int start;   // Start coordinate (x for horizontal, y for vertical)
    int end;     // End coordinate (x for horizontal, y for vertical)
    int height;  // Height (y for horizontal, x for vertical)
    
    ContourSegment(int start, int end, int height)
        : start(start), end(end), height(height) {}
        
    bool operator<(const ContourSegment& other) const {
        return start < other.start;
    }
};

/**
 * Contour class representing the skyline profile of the placement
 * Supports efficient queries and updates in O(log n) time
 */
class Contour {
private:
    // Using a balanced tree (map) to store the contour segments for efficient operations
    std::map<int, int> contourMap;  // Maps coordinate to height/value
    
public:

    Contour();
    
    /**
     * Copy constructor
     */
    Contour(const Contour& other);
    ~Contour();
    void clear();
    
    /**
     * Adds a segment to the contour
     * 
     * @param start Start coordinate of the segment
     * @param end End coordinate of the segment
     * @param height Height of the segment
     */
    void addSegment(int start, int end, int height);
    
    /**
     * Gets the height of the contour at a specific range
     * 
     * @param start Start coordinate of the range
     * @param end End coordinate of the range
     * @return Maximum height in the range
     */
    int getHeight(int start, int end) const;
    
    /**
     * Gets all contour segments
     * 
     * @return Vector of contour segments
     */
    std::vector<ContourSegment> getSegments() const;
    
    /**
     * Merges this contour with another contour
     * 
     * @param other Other contour to merge with
     */
    void merge(const Contour& other);
    
    /**
     * Gets the maximum coordinate value in the contour
     * 
     * @return Maximum coordinate value
     */
    int getMaxCoordinate() const;
    
    /**
     * Gets the maximum height value in the contour
     * 
     * @return Maximum height value
     */
    int getMaxHeight() const;
    
    /**
     * Checks if the contour is empty
     * 
     * @return True if the contour is empty, false otherwise
     */
    bool isEmpty() const;
    
    /**
     * Gets the contour map for direct access
     * (for advanced operations)
     * 
     * @return Reference to the internal contour map
     */
    const std::map<int, int>& getContourMap() const;
};