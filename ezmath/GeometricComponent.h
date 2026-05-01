#ifndef GEOMETRICCOMPONENT_H_
#define GEOMETRICCOMPONENT_H_


namespace GeometrySystem
{
	class GeometricComponent
	{
	public:
		GeometricComponent();
		GeometricComponent(const GeometricComponent& oComponent);
		virtual ~GeometricComponent();
		GeometricComponent& operator=(const GeometricComponent& oComponent);
		virtual void Reset();
		void SetID(const unsigned int& iID);
		unsigned int GetID() const;
		void SetCategory(const int& iCategory);
		int GetCategory() const;

	private:

	protected:
		virtual void Initialize();
		unsigned int m_iID;
		int m_iCategory;
	};
}

#endif


