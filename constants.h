using TypeForM = uint16_t;
const int LIM_BITM = 8 * sizeof(TypeForM);
enum LogicType {ERROR, INPUT, CONSTANT, NOT, AND, OR, XOR, NAND, NOR, XNOR, ANDNY, ANDYN, ORNY, ORYN, MUX };
const uint32_t THRESHOLD = 3;
const uint32_t LENGTH = 15;
const int minimum_lambda = 128;
