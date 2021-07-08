#include <SofaCaribou/Visitor/MultiVecEqualVisitor.h>

namespace SofaCaribou::visitor {
using namespace sofa::core;

auto MultiVecEqualVisitor::fwdMechanicalState(VisitorContext*, behavior::BaseMechanicalState* mm) -> Result {
    if (!p_only_mapped) {
        copy(mm);
    }
    return RESULT_CONTINUE;
}

auto MultiVecEqualVisitor::fwdMappedMechanicalState(VisitorContext*, behavior::BaseMechanicalState* mm) -> Result {
    if (p_mapped || p_only_mapped) {
        copy(mm);
    }

    return RESULT_CONTINUE;
}

auto MultiVecEqualVisitor::getInfos() const -> std::string {
    return "a := b    (with a = '" + p_a.getName() + "', and b = '" + p_b.getName() + "')";
}

template <typename VectorA, typename VectorB>
void internal_copy(VectorA & a, const VectorB & b) {
    if (b.size() > a.size()) {
        a.resize(b.size());
    }

    for (unsigned int i=0; i<b.size(); i++) {
        a[i] = b[i];
    }
}

template <typename DataTypes>
void internal_copy_id(TVecId<V_ALL, V_WRITE> a, TVecId<V_ALL, V_READ> b, State<DataTypes> * state) {
    using namespace sofa::helper;
    using namespace sofa::core;
    using namespace sofa::core::objectmodel;
    using DataVecCoord = Data< typename DataTypes::VecCoord>;
    using DataVecDeriv = Data< typename DataTypes::VecDeriv>;

    if (state == nullptr) {
        throw std::runtime_error("Could not copy b to a.");
    }
    if (a.type == V_COORD) {
        auto a_vector = WriteOnlyAccessor<DataVecCoord>(*state->write(VecCoordId(a)));
        if (b.type == V_COORD) {
            const auto b_vector = ReadAccessor<DataVecCoord>(*state->read(ConstVecCoordId(b)));
            internal_copy(a_vector, b_vector);
        } else {
            // b.type = V_DERIV
            const auto b_vector = ReadAccessor<DataVecDeriv>(*state->read(ConstVecDerivId(b)));
            internal_copy(a_vector, b_vector);
        }
    } else {
        // a.type = V_DERIV
        auto a_vector = WriteOnlyAccessor<DataVecDeriv>(*state->write(VecDerivId(a)));
        if (b.type == V_COORD) {
            const auto b_vector = ReadAccessor<DataVecCoord>(*state->read(ConstVecCoordId(b)));
            internal_copy(a_vector, b_vector);
        } else {
            // b.type = V_DERIV
            const auto b_vector = ReadAccessor<DataVecDeriv>(*state->read(ConstVecDerivId(b)));
            internal_copy(a_vector, b_vector);
        }
    }
}

void MultiVecEqualVisitor::copy(behavior::BaseMechanicalState *state) const {
    using namespace sofa::defaulttype;
    const auto a = p_a.getId(state);
    const auto b = p_b.getId(state);

    if (state->getTemplateName() == Vec1Types::Name()) {
        internal_copy_id(a, b, dynamic_cast<State<Vec1Types>*>(state));
    } else if (state->getTemplateName() == Vec2Types::Name()) {
        internal_copy_id(a, b, dynamic_cast<State<Vec2Types>*>(state));
    } else if (state->getTemplateName() == Vec3Types::Name()) {
        internal_copy_id(a, b, dynamic_cast<State<Vec3Types>*>(state));
    } else {
        throw std::runtime_error("Copy from a State<" + state->getTemplateName() + "> is not yet implemented.");
    }
}

} // namespace SofaCaribou::visitor