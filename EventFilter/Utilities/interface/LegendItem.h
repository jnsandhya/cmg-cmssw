/*
 * LegendItem.h
 *
 *  Created on: Sep 24, 2012
 *      Author: aspataru
 */

#ifndef LEGENDITEM_H_
#define LEGENDITEM_H_

#include <string>

namespace jsoncollector {
class LegendItem {

public:

	LegendItem(std::string name, std::string operation) :
	name_(name), operation_(operation) {}
	virtual ~LegendItem() {}

	std::string const& getName() const {
		return name_;
	}
	std::string const& getOperation() const {
		return operation_;
	}

private:
	std::string name_, operation_;
};
}

#endif /* LEGENDITEM_H_ */
