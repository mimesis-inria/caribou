#include <sofa/simulation/MechanicalVisitor.h>

namespace SofaCaribou::visitor {
/**
 * This visitor is used to compute a = b, where a and b are multi-vectors (MultiVecDeriv or MultiVecCoord).
 * The difference between this visitor and the sofa::simulation::MechanicalVOpVisitor is that here,
 * the assignment from a MutliVecCoord to a MultiVecDeriv is allowed.
 */
class MultiVecEqualVisitor : public sofa::simulation::BaseMechanicalVisitor {
public:

    /**
     * Assign the values of the multi-vector b to the values of a (ie, a := b), no mather
     * the type of a and b (VecDeriv or VecCoord).
     * @param params The current execution parameters
     * @param a The multi-vector a
     * @param b The multi-vector b
     * @param mapped If true, the values of the vectors from mapped mechanical objects are also copied.
     * @param only_mapped If try, only the values of the vectors from mapped mechanical objects are copied.
     */
    MultiVecEqualVisitor(const sofa::core::ExecParams *params,
                         sofa::core::TMultiVecId<sofa::core::V_ALL, sofa::core::V_WRITE> a,
                         sofa::core::TMultiVecId<sofa::core::V_ALL, sofa::core::V_READ> b,
                         bool mapped = false,
                         bool only_mapped = false)
        : BaseMechanicalVisitor(params), p_a(a), p_b(b), p_mapped(mapped), p_only_mapped(only_mapped) {
    }

    // If mapped or only_mapped is ste, this visitor must go through all mechanical mappings, even if isMechanical flag is disabled
    bool stopAtMechanicalMapping(sofa::simulation::Node * /*node*/, sofa::core::BaseMapping *map) override {
        if (p_mapped || p_only_mapped)
            return false;
        else
            return !map->areForcesMapped();
    }

    Result fwdMechanicalState(VisitorContext *ctx, sofa::core::behavior::BaseMechanicalState *mm) override;

    Result fwdMappedMechanicalState(VisitorContext *ctx, sofa::core::behavior::BaseMechanicalState *mm) override;

    void copy(sofa::core::behavior::BaseMechanicalState * state) const ;

    const char *getClassName() const override { return "MultiVecEqualVisitor"; }

    std::string getInfos() const override;

private:
    sofa::core::TMultiVecId<sofa::core::V_ALL, sofa::core::V_WRITE> p_a;
    sofa::core::TMultiVecId<sofa::core::V_ALL, sofa::core::V_READ> p_b;
    bool p_mapped;
    bool p_only_mapped;
};

} // namespace SofaCaribou::visitor
