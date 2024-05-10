#pragma once

#include "gmrf_denoising.hpp"
#include "mylib.hpp"

class processorBase
{
public:
	processorBase() {
		this->prepareToProceess();
	}
	virtual ~processorBase() = 0 {};

	// =================================================================
	virtual void prepareToProceess() = 0 {};
	virtual void process() = 0 {};

	// =================================================================
private:
};

namespace SVGMRF
{
	class processor : public processorBase
	{
		// =================================================================
	public:
		processor() {};
		~processor() override {};

		// =================================================================
		void prepareToProceess() override {};
		void process() override;

		// =================================================================
	private:
	};
}

namespace GMRFUpScale
{
	class processor : public processorBase
	{
		// =================================================================
	public:
		processor() {};
		~processor() override {};

		// =================================================================
		void prepareToProceess() override {};
		void process() override;

		// =================================================================
	private:
	};
}

namespace DVGMRF
{
	class processor : public processorBase
	{
		// =================================================================
	public:
		processor() {};
		~processor() override {};

		// =================================================================
		void prepareToProceess() override {};
		void process() override;

		// =================================================================
	private:
	};
}

namespace SVODGMRF
{
	class processor : public processorBase
	{
		// =================================================================
	public:
		processor() {};
		~processor() override {};

		// =================================================================
		void prepareToProceess() override {};
		void process() override;

		// =================================================================
	private:
	};
}

namespace DVODGMRF {
	class processor : public processorBase
	{
		// =================================================================
	public:
		processor() {};
		~processor() override {};

		// =================================================================
		void prepareToProceess() override {};
		void process() override;

		// =================================================================
	private:
	};
}

namespace DVSVCOMP
{
	class processor : processorBase
	{
		// =================================================================
	public:
		processor() {};
		~processor() override {};

		// =================================================================
		void prepareToProceess() override {};
		void process() override;

		// =================================================================
	private:
	};
}

namespace DVSVODCOMP
{
	class processor : processorBase
	{
		// =================================================================
	public:
		processor() {};
		~processor() override {};

		// =================================================================
		void prepareToProceess() override {};
		void process() override;

		// =================================================================
	private:
	};
}

namespace IVHGMRF {
	class processor : processorBase
	{
		// =================================================================
	public:
		processor() {};
		~processor() override {};

		// =================================================================
		void prepareToProceess() override {};
		void process() override;

		// =================================================================
	private:
	};
}

namespace HGMRFCOMP {
	class processor : processorBase
	{
		// =================================================================
	public:
		processor() {};
		~processor() override {};

		// =================================================================
		void prepareToProceess() override {};
		void process() override;

		// =================================================================
	private:
	};
}

namespace DVHGMRF
{
	class processor : processorBase
	{
		// =================================================================
	public:
		processor() {};
		~processor() override {};

		// =================================================================
		void prepareToProceess() override {};
		void process() override;

		// =================================================================
	private:
	};
}

namespace GMRFCOMP
{
	class processor : processorBase
	{
		// =================================================================
	public:
		processor() {};
		~processor() override {};

		// =================================================================
		void prepareToProceess() override {};
		void process() override;

		// =================================================================
	private:
	};
}