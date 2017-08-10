#include "method.h"

Method* Method::create(int methodCode) {
    switch (methodCode) {
        case METHOD_CODE_HEAT:
            return new MethodHeat();
        case METHOD_CODE_GAS:
            return new MethodGas();
    }
}

