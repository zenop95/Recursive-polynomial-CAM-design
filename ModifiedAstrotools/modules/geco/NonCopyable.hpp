#ifndef geco_NONCOPYABLE_H
#define geco_NONCOPYABLE_H

namespace geco
{
    class NonCopyable
    {
        protected:
            NonCopyable() = default;
            ~NonCopyable() = default;
        private:
            NonCopyable(const NonCopyable&) = delete;
            NonCopyable& operator=(const NonCopyable&) = delete;
    };
}

#endif // geco_NONCOPYABLE_H